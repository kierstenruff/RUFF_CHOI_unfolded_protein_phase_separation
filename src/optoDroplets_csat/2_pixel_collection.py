import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io
import functools
#import napari


from GEN_Utils.FileHandling import df_to_excel
from loguru import logger

logger.info('Import OK')

# define location parameters
image_folder = f'opto_analysis/python_results/initial_cleanup/'
mask_folder = f'opto_analysis/python_results/cellpose_masking/'
output_folder = f'opto_analysis/python_results/pixel_collection/'
if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# define function to collect pixel location and intensity for a given mask, image combination
def pixel_collector(image_array, mask, mask_type=None, visualise=False):
    pixel_array = np.where(mask == 1, image_array, np.nan)
    # plt.imshow(pixel_array)
    coords = pd.DataFrame(pixel_array).unstack().reset_index().dropna()
    coords.columns = ['x', 'y', 'intensity']

    if mask_type != None:
        coords['mask_type'] = mask_type

    if visualise:
        # test visualisation, compare to plt.show
        fig, ax = plt.subplots(figsize=(20, 20))
        sns.scatterplot(coords['x'], coords['y'], hue=coords['intensity'], palette='magma_r', size=0.5,linewidth=0, alpha = 0.7)
        plt.ylim(1024, 0)
        plt.xlim(0, 1024)

    return coords


# read in images, including transpose to put in x, y, z format
file_list = [filename for filename in os.listdir(image_folder) if '.tif' in filename]
images = {filename.replace('.lif - ', '_').replace('.tif', ''): skimage.io.imread(f'{image_folder}{filename}').transpose(1, 2, 0) for filename in file_list}

# read in masks - remember that stack format is [whole_cell, nucleus, cytoplasm], 
# fails try/except if no masks found therefore skip that image
masks = {}
for image_name, img in images.items():
    logger.info(f'Processing {image_name}')
    try:
       masks[image_name] = {cell_number.replace('.npy', ''): np.load(f'{mask_folder}{image_name}/{cell_number}') for cell_number in os.listdir(f'{mask_folder}{image_name}')}
       logger.info(f'Masks loaded for {len(masks[image_name].keys())} cells')
    except:
        logger.info(f'{image_name} not processed as no mask found')


# # Example napari visualisation to test matching mask array to image
# image_test_name = 'Hsp70_I88G-NLS_1'
# with napari.gui_qt():
#     viewer = napari.Viewer()
#     viewer.add_image(images[image_test_name][:, :, 0], name='raw_image')
#     for cell_number, array_stack in masks[image_test_name].items():
#         array_stack
#         viewer.add_labels(array_stack[:, :, 0], name=f'{cell_number}')

# collect pixel information
pixel_information = {}
for image_name, img in images.items():
    logger.info(f'Processing {image_name}')
    try:
        channels = []
        for image_channel in range(img.shape[2]):
            logger.info(f'Processing channel {image_channel} pixels')
            # collect only one channel
            image_array = img[:, :, image_channel]
            # for each cell, collect cell, nuc, cyto pixels
            cells = []
            for cell_number, mask_array in masks[image_name].items():
                cell_number
                # collect pixels and coords
                cell_pixels = pd.concat([pixel_collector(image_array, mask_array[:, :, i], visualise=False, mask_type=mask_type) for i, mask_type in enumerate(['cell', 'nucleus', 'cytoplasm'])])
                # add identifiers
                cell_pixels['cell_number'] = cell_number
                cell_pixels.rename(columns={'intensity': f'channel_{image_channel}'}, inplace=True)
                cells.append(cell_pixels)
            channels.append(pd.concat(cells))
        channels = functools.reduce(lambda left, right: pd.merge(left, right, on=['x', 'y', 'cell_number', 'mask_type'], how='outer'), channels)
        channels['image_name'] = image_name
            
        pixel_information[image_name] = channels
    except:
        logger.info(f'{image_name} not processed as no mask found')
logger.info('Completed pixel collection')

# save to excel (although this will likely be unopenable if more than a few images) and csv
# df_to_excel(output_path=f'{output_folder}pixel_information.xlsx', sheetnames=list(pixel_information.keys()), data_frames=list(pixel_information.values()))
saved = [df.to_csv(f'{output_folder}{image_name}.csv') for image_name, df in pixel_information.items()]
