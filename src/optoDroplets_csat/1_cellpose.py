import os
import numpy as np
import matplotlib.pyplot as plt
import skimage.io
from cellpose import models
from cellpose import plot
import collections
#import napari
from loguru import logger

input_folder = f'opto_analysis/python_results/initial_cleanup/'
output_folder = f'opto_analysis/python_results/cellpose_masking/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

def apply_cellpose(images, image_type='cyto', channels=[0, 0], diameter=None, flow_threshold=0.4, cellprob_threshold=0.0):
    """Apply model to list of images. Returns masks, flows, styles, diams.
    - model type is 'cyto' or 'nuclei'
    - define CHANNELS to run segementation on (grayscale=0, R=1, G=2, B=3) where channels = [cytoplasm, nucleus]. If NUCLEUS channel does not exist, set the second channel to 0
    """
    model = models.Cellpose(model_type=image_type)
    masks, flows, styles, diams = model.eval(
        images, diameter=diameter, channels=channels, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold)
    return masks, flows, styles, diams
    
def visualise_cell_pose(images, masks, flows, channels=[0,0]):
    """Display cellpose results for each image
    """
    for image_number, image in enumerate(images):
        maski = masks[image_number]
        flowi = flows[image_number][0]

        fig = plt.figure(figsize=(12,5))
        plot.show_segmentation(fig, image, maski, flowi, channels=channels)
        plt.tight_layout()
        plt.show()

def edge_filter(mask):
    """Collect boundary pixel values for all edges, return unique values
    which correspond to cells that are touching/over the edge boundaries""" 
    size = mask.shape[0]
    edge_1_cells = set(list(mask[0:1, 0:size].flatten()))
    edge_2_cells = set(list(mask[0:size, (size-1):size].flatten()))
    edge_3_cells = set(list(mask[(size-1):size, 0:size].flatten()))
    edge_4_cells = set(list(mask[0:size, 0:1].flatten()))
    return edge_1_cells | edge_2_cells | edge_3_cells | edge_4_cells

def size_filter(mask, lower_size=1500, upper_size=10000):
    """Collect cells that are outside the cell size bounds as those to
    be excluded""" 
    cell_size = dict(collections.Counter(mask.flatten()))
    bs_cells = [
        cell_number
        for cell_number, cell_size in cell_size.items()
        if cell_size < lower_size or cell_size > upper_size
    ]

    return set(bs_cells)

# --------------------------------------Initialise file list--------------------------------------

file_list = []
for filename in os.listdir(input_folder):
    if filename.endswith('.tif'):
        file_list.append(filename)
    else:
        continue 
# file_list = [filename for filename in os.listdir(input_folder)]
#with napari.gui_qt():
 #   viewer = napari.view_image(skimage.io.imread(f'{input_folder}{file_list[0]}'))

# reading in all channels for each image, and transposing to correct dimension of array
imgs = [skimage.io.imread(f'{input_folder}{filename}').transpose(1, 2, 0) for filename in file_list]

# clean filenames
img_names = [filename.replace('.tif', '') for filename in file_list]

# -----------------------Complete cellpose with cytoplasm channel---------------------------------
#Channel 0: Hoerscht
#Channel 1: BF
#Channel 2: Antibody - this looks ok, and certainly better than brightfield. Need to be careful for those that are clumped together - may need to consider Bf

# collecting only channel 0's for cytoplasm
cytoplasm_images = [image[:, :, 0] for image in imgs]
plt.imshow(cytoplasm_images[0])

# Apply cellpose then visualise
masks, flows, styles, diams = apply_cellpose(cytoplasm_images, image_type='cyto', diameter=80, flow_threshold=0.8)
visualise_cell_pose(cytoplasm_images, masks, flows, channels=[0, 0])

# filter identified cell number according to boundary and cell size (reasonable range ~1500, 10000)
# using set logic here - exclude cells will be unique set of cells
exclude_cells = {}
for image_number, mask in enumerate(masks):
    exclude_cells[image_number]  = list(edge_filter(mask) | size_filter(mask, lower_size=20, upper_size=100000))

# remove cells from mask according to exclude cells:
filtered_masks = {}
for image_number, exclude_vals in exclude_cells.items():
    mask = masks[image_number].copy()
    plt.imshow(mask)
    cell_mask = np.where(~np.isin(mask, exclude_vals), mask, np.nan)
    plt.imshow(cell_mask)
    filtered_masks[image_number] = cell_mask

# -----------------------Complete cellpose with nuclei channel---------------------------------
nuclei_images = [image[:, :, 8] for image in imgs]
# plt.imshow(cytoplasm_images[0]+nuclei_images[0])

# Apply cellpose then visualise
nuc_masks, nuc_flows, nuc_styles, nuc_diams = apply_cellpose(nuclei_images, diameter=50, image_type='nuclei')
visualise_cell_pose(nuclei_images, nuc_masks, nuc_flows, channels=[0, 0])

# save associated cell mask arrays
np.save(f'{output_folder}cellpose_masks.npy', masks)
np.save(f'{output_folder}cellpose_nuclei.npy', nuc_masks)
np.save(f'{output_folder}cellpose_filtered_masks.npy', np.dstack(filtered_masks.values()))
# ------------------------Match nuclei to cells-----------------------------------------------
#filtered_masks = np.load(f'{output_folder}cellpose_filtered_masks.npy')
nuc_masks = np.load(f'{output_folder}cellpose_nuclei.npy')
#filtered_masks = {image_number: filtered_masks[:, :, image_number] for image_number in range(filtered_masks.shape[2])}
final_masks = {}
for image_number, mask in filtered_masks.items():
    logger.info(f'processing image number {image_number}')
    cyto_mask = mask.copy()
    nuc_mask = nuc_masks[image_number].copy()
    #plt.imshow(cyto_mask+nuc_mask)
    for cell_number in np.unique(cyto_mask[~np.isnan(cyto_mask)]):
        logger.info(f'processing cell number {cell_number}')
        whole_cell = np.where(cyto_mask == cell_number, 1, 0)
        nucleus = np.where(whole_cell == 1, nuc_mask, 0)
        # to filter nuclei for all within cell, check count nuc_id before + after filtering
        nuc_numbers = [val for val in list(np.unique(nucleus)) if val != 0]
        if not len(nuc_numbers) == 1:
            continue
        nucID = nuc_numbers[0]
        nucID_count_CELL = nucleus.flatten().tolist().count(nucID)
        nucID_count_TOTAL = nuc_mask.flatten().tolist().count(nucID)
        if not nucID_count_CELL == nucID_count_TOTAL:
            continue
        cytoplasm = np.where(nucleus == 0, whole_cell, 0)
        nucleus = np.where(nucleus != 0, 1, 0)
        final_masks[(image_number, cell_number)] = np.dstack([whole_cell, nucleus, cytoplasm])

# Save individual filtered cell masks
for (image_number, cell_number), array_stack in final_masks.items():
    image_name = img_names[image_number]
    # create folder for each image output
    if not os.path.exists(f'{output_folder}{image_name}/'):
        os.mkdir(f'{output_folder}{image_name}/')
    # save associated cell mask arrays
    np.save(f'{output_folder}{image_name}/cell_{int(cell_number)}.npy', array_stack)
for image_name, image in zip(img_names, imgs):
    # save original image
    np.save(f'{output_folder}{image_name}/{image_name}.npy', image)
# optional - edit nuclei masks?
# nuc_masks = segmentation_editor(images, masks)

