#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Processing raw images to remove / correct artifacts
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

                                           
###############################################################################
## Command Line Interface
###############################################################################
from MultiplexImaging_utils import *


argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input-dir", help="Input directory containing TIFF files.", type=dir_path)
argParser.add_argument("-o", "--output-dir", help="Output Directory", type=dir_path)
argParser.add_argument("-s", "--samples", help="Sample name to look for.", type=str)
argParser.add_argument("-m", "--marker-metadata", help="Marker annotation containing QC directives.", type=str)
argParser.add_argument("-a", "--annotation-dir", help="Annotation directory containing potential regions to remove.", type=dir_path)

args = argParser.parse_args()

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
from skimage.segmentation import find_boundaries
from skimage.morphology import remove_small_objects, h_maxima, remove_small_holes
from matplotlib import pyplot as plt
import numpy as np
import time
from skimage.segmentation import expand_labels
import math


print()
print("Running QC Filtering...")
print()
print("######################################################################################")     
print("Input directory = %s" % args.input_dir)
print("Output directory = %s" % args.output_dir)
samples = [str(s.strip()) for s in args.samples.split(",")]
samples.pop()
print("Samples = %s" % samples)
print("Marker metadata = %s" % args.marker_metadata)
marker_metadata = pd.read_csv(os.path.join(args.marker_metadata))
print("Annotation directory = %s" % args.annotation_dir)


markers = np.array(marker_metadata.Marker[marker_metadata.PassOverallQuality == True])
print(markers)

spot_artifact_markers = np.array(marker_metadata.Marker[marker_metadata.RemoveSpots == True])
print("Markers for which to detect and remove bright artifact spots = %s" % spot_artifact_markers)


marker_metadata.RemovePixelBelow_Samples = marker_metadata.RemovePixelBelow_Samples.apply(str) 

spills_artifact_dict = {}
for marker in np.array(marker_metadata.Marker[(marker_metadata.RemovePixelBelow_Samples != "nan") == True]):
    print(marker)
    df = marker_metadata[marker_metadata.Marker == marker]
    samples_to_filter = np.array(df.RemovePixelBelow_Samples.values[0].split(';'))
    thresholds = np.array(df.RemovePixelBelow_Thresholds.values[0].split(';'))
    spills_artifact_dict[marker] = {}
    for i in range(len(samples_to_filter)):
        threshold = thresholds[i]
        sample = samples_to_filter[i]
        spills_artifact_dict[marker].update({sample : int(threshold)})
    print(spills_artifact_dict[marker])


output_dir = os.path.join(args.output_dir)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
tiff_dir = os.path.join(output_dir, "input")
if not os.path.exists(tiff_dir):
    os.mkdir(tiff_dir)

start_time = time.time()

print(spills_artifact_dict)
###############################################################################
## QC Filtering 
###############################################################################

for sample in samples:
    print("Running sample: " + sample)
    
    print("######################################################################################")    
    print()
    print("Copying and cropping images to the right dimensions...")
    initialize_input(args.input_dir, output_dir, sample, markers)


    napari_shape_to_remove = os.path.join(args.annotation_dir, sample + "_remove_mask.csv")
    if os.path.exists(napari_shape_to_remove):
        print("######################################################################################")    
        print()
        print("Removing artifact regions: manually identified noisy regions...")
    
        image = tif.imread(os.path.join(tiff_dir, sample, markers[0] + ".tiff" ))
        napari_polygon_to_mask(sample, image.shape, napari_shape_to_remove, tiff_dir, mask_suffix = "_remove_mask.tiff")
        mask = os.path.join(tiff_dir, sample, sample + "_remove_mask.tiff")
        remove_artifacts_from_single_mask(mask, tiff_dir, [sample], markers)
    
    if len(spot_artifact_markers) > 0:
        print("######################################################################################")    
        print()
        print("Removing artifact regions: automatic detection of bright spots...")
        detect_bright_artifacts(tiff_dir, [sample], markers = spot_artifact_markers, mask_suffix = "_spots.tiff") 
        remove_artifacts_from_multiple_masks(tiff_dir, [sample], markers = spot_artifact_markers, mask_suffix = "_spots.tiff")
        
        
    if len(spills_artifact_dict) > 0:
        print("######################################################################################")
        print()
        print("Removing artifact regions: removing spills based on marker metadata...")
        spilled_markers = []
        spilled_thresholds = []
        for marker in spills_artifact_dict:
            print(marker)
            if (sample in spills_artifact_dict[marker].keys()):
                print(spills_artifact_dict[marker])
                spilled_markers.append(marker)
                spilled_thresholds.append(spills_artifact_dict[marker][sample])
            if  ("All" in spills_artifact_dict[marker].keys()):
                spilled_markers.append(marker)
                spilled_thresholds.append(spills_artifact_dict[marker]["All"])
        print(spilled_markers)
        print(spilled_thresholds)
        remove_artifacts_below_threshold(spilled_thresholds, tiff_dir, sample, spilled_markers)
    

    print("######################################################################################")
    print()


print('Finished Running QC Filtering in %s seconds', round(time.time() - start_time, 2))

