import os
import re
from shutil import copyfile


#logger = logger_config(__name__)
#logger.info('Import ok')

def jarvis(input_path, output_path):

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for filename in os.listdir(input_path):
        if filename.endswith('.tif'):
            filename
            pathname = os.path.join(input_path, filename)
            if os.path.isdir(pathname):
                print (f"Directory not processed {filename}")
            else:
                details = re.split('.lif -', filename)
                print(details)
                new_name = details[0]+'_'+details[1]
                copyfile(input_path+filename, output_path+new_name)
        else:
            continue
                
jarvis(input_path= f'opto_analysis/data/', output_path = f'opto_analysis/python_results/initial_cleanup/')
