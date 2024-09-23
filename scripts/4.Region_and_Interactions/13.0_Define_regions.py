#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Simple region detection using custom functions
# Finding :
# - Endothelial vessel regions based on Actin
# - Lymphatic vessel regions based on Podoplanin
# - Epidermis regions based on Cytokeratin
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

###############################################################################
# # Command Line Interface
# ##############################################################################
from MultiplexImaging_utils import *

argParser = argparse.ArgumentParser()
argParser.add_argument("-o", "--output-dir", help="Output Directory", type=dir_path)

args = argParser.parse_args()

print()
print("Running Region detection...")
print()
print("######################################################################################")     
print("Output directory = %s" % args.output_dir)



output_dir = os.path.join(args.output_dir, "Regions")
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
tiff_dir = os.path.join(args.output_dir, "input")
samples = os.listdir(tiff_dir)
segmentation_dir = os.path.join(output_dir, "segmentation")

marker_file  = "annotation/marker_metadata.csv"
marker_df = pd.read_csv(os.path.join(marker_file))

###############################################################################
# # Paths & Loading
# ############################################################################## 
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
from skimage.segmentation import find_boundaries
from skimage.morphology import remove_small_objects, h_maxima, remove_small_holes
from matplotlib import pyplot as plt
import numpy as np
import time
from skimage.segmentation import expand_labels
import pandas as pd
import numpy as np
from tifffile import imread
import random
import re
import matplotlib.pyplot as plt
from skimage import data
from skimage.filters import threshold_otsu
import numpy as np
import matplotlib.pyplot as plt
from skimage import data, img_as_float
from skimage.segmentation import (morphological_chan_vese,
                                  morphological_geodesic_active_contour,
                                  inverse_gaussian_gradient,
                                  checkerboard_level_set)
from skimage.measure import label, regionprops, regionprops_table
from skimage.draw import ellipse, polygon2mask
from skimage.morphology import reconstruction
from skimage.filters import gaussian
from tifffile import imwrite
start_time = time.time()

###############################################################################
# # Preprocessing
# ##############################################################################

object_sizes = [40000, 5000, 15000]
markers = ["Cytokeratin", "Actin", "Podoplanin"]
add = [0, 1000, 2000]

for i in [0, 1, 2]:
    marker = markers[i]
    for sample in samples:

        contrast_min_marker = marker_df.ContrastRange_min[marker_df.Marker == marker].values[0]
        contrast_max_marker = marker_df.ContrastRange_min[marker_df.Marker == marker].values[0]

        image = imread(os.path.join(tiff_dir, sample, marker + ".tiff"))
        image[image < (contrast_min_marker + 0.25 * contrast_min_marker)] = 0
        image[image < contrast_max_marker] = contrast_max_marker + 0.25 * contrast_max_marker

        quant = np.quantile(image[image > 0], 0.99)
        image[image > quant] = quant
        
        image2 = gaussian(image, sigma=20)
        thresh = threshold_otsu(image2)
        binary = image2 > thresh
        seed = np.copy(binary)
        seed[1:-1, 1:-1] =binary.max()
        mask = binary
        filled = reconstruction(seed, mask, method='erosion')

        filled = filled.astype(int)
        filled = label(filled)
        regions = regionprops(filled)
        print("Regions = " + str(len(regions)))
        for j in range(len(regions)):
            if regions[j].area <  object_sizes[i]:
                filled[filled == regions[j].label] = 0
        
        filled[filled>0]= filled[filled>0] + add[i]
        filled = filled.astype(int)
        
        imwrite(os.path.join(output_dir, sample + "-" + marker + "-mask.tiff"), filled)



for sample in samples:
    print(sample)
    cyto = imread(os.path.join(output_dir, sample + "-Cytokeratin-mask.tiff") )
    actin = imread(os.path.join(output_dir, sample + "-Actin-mask.tiff"))
    podo = imread(os.path.join(output_dir, sample + "-Podoplanin-mask.tiff"))
    combined = podo
    combined[actin > 0] = actin[actin > 0]
    combined[cyto > 0] = cyto[cyto > 0]
    imwrite(os.path.join(output_dir, sample + "-Cytokeratin-mask.tiff"), combined)
    
print('Finished Running Define regions in %s seconds', round(time.time() - start_time, 2))

