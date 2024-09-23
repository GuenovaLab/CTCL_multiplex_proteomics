#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:33:31 2023

@author: localadmin
"""

import os
import pandas as pd
import numpy as np
from tifffile import imread
import random
import re
import cv2 as cv

import napari
from napari.utils import nbscreenshot

os.environ["QT_QPA_PLATFORM"] = "wayland-xcomposite-glx"

# %gui

# Create an empty viewer
viewer = napari.Viewer()

# Change below
output_dir = "../output/manual_annotation_marker2/"
tiff_dir = "../output/input/"
segmentation_dir  = "../output/segmentation/"
marker_file  = "../annotation/marker_metadata.csv"
marker_df = pd.read_csv(os.path.join(marker_file))
marker_unique = marker_df.Marker[(marker_df.Marker != "DAPI") & (marker_df.PassOverallQuality == True)]
marker_unique = marker_unique.values

