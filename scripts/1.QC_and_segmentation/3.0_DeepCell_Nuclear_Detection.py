#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Cell segmentation using Mesmer (Deepcell)
# Nuclear marker = DAPI
# Cytoplasmic marker = combination of Cytokeratin, CD45RO, betaActin, Galectin3, p16INK4a
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

# Scripts adapted from ARK analysis pipeline:
# https://github.com/angelolab/ark-analysis
#
# Greenwald, N.F., Miller, G., Moen, E. et al. Whole-cell segmentation of tissue 
# images with human-level performance using large-scale data annotation and deep learning.
# Nat Biotechnol 40, 555–565 (2022). https://doi.org/10.1038/s41587-021-01094-0
#
# Liu, C.C., Greenwald, N.F., Kong, A. et al. Robust phenotyping of highly multiplexed
# tissue imaging data using pixel-level clustering. Nat Commun 14, 4618 (2023).
# https://doi.org/10.1038/s41467-023-40068-5 

###############################################################################
## Command Line Interface
###############################################################################
from MultiplexImaging_utils import *
#os.environ['LD_LIBRARY_PATH'] = "/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cuda_runtime/lib/"

argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input-dir", help="Input directory containing TIFF files.", type=dir_path)
argParser.add_argument("-s", "--sample-name", help="Sample name to look for.", type=str)
argParser.add_argument("-o", "--output-dir", help="Output Directory", type=dir_path)
argParser.add_argument("-r", "--resolution", help="Image resolution in um/pixel", type=float, default=0.17)
argParser.add_argument("-n", "--nucleus-channels", help="Channel for nucleus channel, separated by commas. Defaults to DAPI.", type=str, default="DAPI")
argParser.add_argument("-c", "--cytoplasm-channels", help="Input Cytoplasm channel, separated by comma. Defaults to Cytokeratin,CD45RO,betaActin,Galectin3,p16INK4a.", type=str, default="Cytokeratin,CD45RO,betaActin,Galectin3,p16INK4a")
argParser.add_argument("-C", "--postprocess-args-cytoplasm", metavar="KEY=VALUE", nargs='+', help="Arguments passed to postprocess_kwargs_whole_cell. In the form param=value. \
                       Parameters include interior_threshold, maxima_threshold...")
argParser.add_argument("-N", "--postprocess-args-nucleus", metavar="KEY=VALUE", nargs='+', help="Arguments passed to postprocess_kwargs_nuclear. In the form param=value. \
                       Parameters include pixel_expansion...")
argParser.add_argument("-m", "--marker-metadata", help="Marker annotation.", type=str)
argParser.add_argument("-t", "--minimum-cell-area", help="Minimum cell area in um. Default to 78.53982um² (10um diameter).", type=float, default = 78.53982)

args = argParser.parse_args()

print()
print("Running Cell Segmentation with DeepCell...")
print()
print("######################################################################################")     
print("Input directory = %s" % args.input_dir)
print("Output directory = %s" % args.output_dir)
print("Sample Name = %s" % args.sample_name)
print("Resolution = %s" % args.resolution)
print("Nuclear channel = %s" % args.nucleus_channels)
print("Cytoplasmic channel = %s" % args.cytoplasm_channels)
print("Postprocess nucleus arguments = %s" % args.postprocess_args_nucleus)
print("Postprocess cytoplasm arguments = %s" % args.postprocess_args_cytoplasm)
print("Marker metadata = %s" % args.marker_metadata)
print("Minimum cell area = %s" % args.minimum_cell_area)
marker_metadata = pd.read_csv(os.path.join(args.marker_metadata))


output_dir = os.path.join(args.output_dir)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
tiff_dir = os.path.join(output_dir, "input")
if not os.path.exists(tiff_dir):
    os.mkdir(tiff_dir)
segmentation_dir = os.path.join(output_dir, "segmentation")
if not os.path.exists(segmentation_dir):
    os.mkdir(segmentation_dir)
cell_table_dir = os.path.join(output_dir, "cell_table")
if not os.path.exists(cell_table_dir):
    os.mkdir(cell_table_dir)
     
markers = np.array(marker_metadata.Marker[marker_metadata.PassOverallQuality == True])
print(markers)

nucleus_channels = [str(s.strip()) for s in args.nucleus_channels.split(",")]
cytoplasm_channels = [str(s.strip()) for s in args.cytoplasm_channels.split(",")]

print("nucleus_channels %s"  % nucleus_channels)
print("cytoplasm_channels %s"  % cytoplasm_channels)

nucleus_files = dict()
for i in nucleus_channels:
    nucleus_files[i] = os.path.join(tiff_dir, args.sample_name, i + ".tiff")

cytoplasm_files = dict()
for i in cytoplasm_channels:
    cytoplasm_files[i] = os.path.join(tiff_dir, args.sample_name, i + ".tiff")


print("Nucleus image = %s" % nucleus_files)
print("Cytoplasm images = %s" % cytoplasm_files)


postprocess_kwargs_nuclear = {}
if args.postprocess_args_nucleus is not None:
    postprocess_kwargs_nuclear = parse_vars(args.postprocess_args_nucleus, "float")
    print(postprocess_kwargs_nuclear)

postprocess_kwargs_whole_cell = {}
if args.postprocess_args_cytoplasm is not None:
    postprocess_kwargs_whole_cell = parse_vars(args.postprocess_args_cytoplasm, "float")
    print(postprocess_kwargs_whole_cell)

###############################################################################
## Paths & Loading
############################################################################### 
import sys
from PIL import Image
import tifffile as tif
from alpineer import image_utils, io_utils, misc_utils
from deepcell.utils.plot_utils import create_rgb_image
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import make_outline_overlay
from ark.segmentation import marker_quantification
from ark.utils import (data_utils, deepcell_service_utils, example_dataset,
                       plot_utils)
from ark.segmentation import fiber_segmentation, marker_quantification, regionprops_extraction, segmentation_utils, signal_extraction 
from skimage.segmentation import find_boundaries
from skimage.morphology import remove_small_objects, h_maxima, remove_small_holes
from matplotlib import pyplot as plt
import numpy as np
import time
from skimage.segmentation import expand_labels
import tensorflow as tf
print(tf.config.list_physical_devices("GPU"))

start_time = time.time()

###############################################################################
## Preprocessing
###############################################################################jup
print("Converting to Green and Blue channels...")
nucleus_images = dict()
nucleus = np.array(0)

for i in nucleus_files:
    nucleus_images[i] = Image.open(nucleus_files[i])
    nucleus_images[i] = nucleus_images[i].point(lambda i: i/255).convert("L")
    nucleus = nucleus + np.array(nucleus_images[i])

cytoplasm_images = dict()
cytoplasm = np.array(0)
for i in cytoplasm_files:
    cytoplasm_images[i] = Image.open(cytoplasm_files[i])
    cytoplasm_images[i] = cytoplasm_images[i].point(lambda i: i/255).convert("L")
    cytoplasm = cytoplasm + np.array(cytoplasm_images[i])

rgb_array =  np.zeros((nucleus.shape[0],nucleus.shape[1],3), 'uint8')
rgb_array[..., 1] = nucleus
rgb_array[..., 2] = cytoplasm
image = Image.fromarray(rgb_array)
rgb_array = rgb_array.reshape((1, image.height, image.width, 3))
rgb_array = rgb_array[:, :, :, 1:3]

print("Image format = %s" % image.format, image.size, image.mode)

rgb_images = create_rgb_image(rgb_array, channel_colors=['green', 'blue'])

# plot the data
fig, ax = plt.subplots(1, 3, figsize=(15, 15))
ax[0].imshow(rgb_array[0, ..., 0])
ax[1].imshow(rgb_array[0, ..., 1])
ax[2].imshow(rgb_images[0, ...])

ax[0].set_title('Nuclear channel')
ax[1].set_title('Membrane channel')
ax[2].set_title('Overlay')

plt.savefig(os.path.join(segmentation_dir, args.sample_name + "_channels.pdf"), bbox_inches='tight')

###############################################################################
## Nuclear Segmentation
###############################################################################
app = Mesmer()
print('Running Mesmer Segmentation...')
segmentation_predictions = app.predict(rgb_array,
                                       image_mpp = args.resolution,
                                       compartment = "both", 
                                       postprocess_kwargs_whole_cell = postprocess_kwargs_whole_cell,
                                       postprocess_kwargs_nuclear = postprocess_kwargs_nuclear)

print("Removing cells smaller than 10um...")
segmentation_predictions = segmentation_predictions.copy()

whole_cell_image = segmentation_predictions.copy()[0, ..., 0]
whole_cell_image_before = Image.fromarray(whole_cell_image)
whole_cell_image_before.save(os.path.join(segmentation_dir,  args.sample_name + "_whole_cell_before_filtering" + ".tiff"))

regions = regionprops(segmentation_predictions.copy()[0, ..., 0])
print("Regions = " + str(len(regions)))
for i in range(len(regions)):
    if regions[i].area <  int(args.minimum_cell_area / args.resolution):
        whole_cell_image[whole_cell_image == regions[i].label] = 0

whole_cell_image = expand_labels(whole_cell_image, distance=10)


nuclear_image = segmentation_predictions.copy()[0, ..., 1]
nuclear_image_before = Image.fromarray(nuclear_image)
nuclear_image_before.save(os.path.join(segmentation_dir,  args.sample_name + "_nuclear_before_filtering" + ".tiff"))

regions = regionprops(segmentation_predictions.copy()[0, ..., 1])
print("Regions = " + str(len(regions)))
   
for i in range(len(regions)):
    if regions[i].area <  int(0.75 * args.minimum_cell_area / args.resolution):
        nuclear_image[nuclear_image == regions[i].label] = 0
        

###############################################################################
## Exporting Restults
###############################################################################

print('Saving Segmentation...')
overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=whole_cell_image.reshape(1,whole_cell_image.shape[0],whole_cell_image.shape[1],1))
overlay_image = 255 * overlay_data[0, ...]
overlay_image = Image.fromarray(overlay_image.astype(np.uint8), "RGB")
overlay_image.save(os.path.join(segmentation_dir, args.sample_name + "_overlay_whole_cell" + ".tiff"))
whole_cell_image = Image.fromarray(whole_cell_image)
whole_cell_image.save(os.path.join(segmentation_dir,  args.sample_name + "_whole_cell" + ".tiff"))

overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=nuclear_image.reshape(1,nuclear_image.shape[0],nuclear_image.shape[1],1))
overlay_image = 255 * overlay_data[0, ...]
overlay_image = Image.fromarray(overlay_image.astype(np.uint8), "RGB")
overlay_image.save(os.path.join(segmentation_dir, args.sample_name + "_overlay_nuclear" + ".tiff"))
nuclear_image = Image.fromarray(nuclear_image)
nuclear_image.save(os.path.join(segmentation_dir,  args.sample_name + "_nuclear" + ".tiff"))


# Saving Whole Cell Borders
labels = segmentation_utils.load_utils.load_imgs_from_dir(data_dir=segmentation_dir,
                                           files=[args.sample_name + '_whole_cell.tiff'],
                                           xr_dim_name='compartments',
                                           xr_channel_names=['whole_cell'],
                                           trim_suffix='_whole_cell',
                                           match_substring='_whole_cell')
# generates segmentation borders and labels
labels = labels.loc[args.sample_name, :, :, 'whole_cell'].values
# define borders of cells in mask
contour_mask = find_boundaries(labels, connectivity=1, mode='inner').astype(np.uint8)
contour_mask[contour_mask > 0] = 255
# save the cell border image
save_path_seg_borders = os.path.join(segmentation_dir, f'{args.sample_name}_whole_cell_segmentation_borders.tiff')
image_utils.save_image(save_path_seg_borders, contour_mask)

# Saving Nuclear Borders
labels = segmentation_utils.load_utils.load_imgs_from_dir(data_dir=segmentation_dir,
                                           files=[args.sample_name + '_nuclear.tiff'],
                                           xr_dim_name='compartments',
                                           xr_channel_names=['nuclear'],
                                           trim_suffix='_nuclear',
                                           match_substring='_nuclear')
# generates segmentation borders and labels
labels = labels.loc[args.sample_name, :, :, 'nuclear'].values
# define borders of cells in mask
contour_mask = find_boundaries(labels, connectivity=1, mode='inner').astype(np.uint8)
contour_mask[contour_mask > 0] = 255
# save the cell border image
save_path_seg_borders = os.path.join(segmentation_dir, f'{args.sample_name}_nuclear_segmentation_borders.tiff')
image_utils.save_image(save_path_seg_borders, contour_mask)

print('Extracting Single-Cell Matrix...')
# set to True to add nuclear cell properties to the expression matrix
nuclear_counts = True

# set to True to bypass expensive cell property calculations
# only cell label, size, and centroid will be extracted if True
fast_extraction = False

cell_table_size_normalized, cell_table_arcsinh_transformed = \
    marker_quantification.generate_cell_table(segmentation_dir=segmentation_dir,
                                              tiff_dir=tiff_dir,
                                              img_sub_folder=None,
                                              fovs=[args.sample_name],
                                              batch_size=5,
                                              nuclear_counts=nuclear_counts,
                                              fast_extraction=fast_extraction)

print('Saving Single-Cell Matrix...')
compression = "gzip"
cell_table_size_normalized.to_csv(os.path.join(cell_table_dir, args.sample_name + '_cell_table_size_normalized.csv.gz'),
                                  compression=compression, index=False)
cell_table_arcsinh_transformed.to_csv(os.path.join(cell_table_dir, args.sample_name + '_cell_table_arcsinh_transformed.csv.gz'),
                                      compression=compression, index=False)


print('Finished Running Segmentation in %s seconds', round(time.time() - start_time, 2))

