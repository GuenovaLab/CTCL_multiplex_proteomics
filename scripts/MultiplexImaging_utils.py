#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:32:50 2023

@author: Pacôme Prompsy
"""
import json, os, re, argparse, feather
from typing import List
from skimage.io import imread
from alpineer import image_utils, io_utils, load_utils, misc_utils
import glob
import pandas as pd
import os
import pathlib
import re
from typing import List, Union
from PIL import Image
import tifffile as tif
import feather
import natsort as ns
import numpy as np
from tqdm.notebook import tqdm_notebook as tqdm
from ark.utils import data_utils
from ark import settings
from skimage import segmentation
from skimage.measure import label, regionprops, regionprops_table
from skimage.draw import ellipse, polygon2mask
from skimage.morphology import reconstruction
from skimage.filters import gaussian
###############################################################################
## Multiplexing specfici util functions
###############################################################################

def read_markers(marker_file):
    with open(marker_file, "r") as fp:
         lines = fp.readlines()
    markers = []
    for l in lines:
        markers.append(l.replace("\n", ""))
    return(markers)
    
def qiSettings_to_markers(qiSettings_file):
    """
    Create input links correctly named and create output directory structure.
    Returns the new input directory containing links to TIFF files.
    """
    f = open(qiSettings_file)
    data = json.load(f)
    markers = []
    for item in data['data']:
        markers.append(item['comparisonData']['dye'].replace(" ", ""))
    return markers

def initialize_input(input_dir, output_dir, sample_name, markers):
    """
    Create input links correctly named and create output directory structure.
    Returns the new input directory containing links to TIFF files.
    """
    
    tiffdir = os.path.join(output_dir, "input", sample_name)
    if not os.path.exists(tiffdir):
        os.makedirs(tiffdir)
            
    TIFFs = os.listdir(os.path.join(input_dir, sample_name))
    links = dict.fromkeys(markers, "")

    for file in TIFFs:
        for marker in markers:
            if marker == "DAPI" or marker == "TRBC1":
                regex = re.compile(".*-"  + marker + ".*.tif.*")
            else:
                regex = re.compile(".*-"  + marker + "_.*.tif.*")
            
            if regex.match(file):
                links[marker] = file

    print("Found this files:")
    print(links)
    
    widthsX = dict()
    heigthY = dict()
    for marker in links:
        url = os.path.join(input_dir, sample_name, links[marker])
        if os.path.isfile(url):    
            img = Image.open(url)
            widthsX[marker] = int(img.size[0])
            heigthY [marker] = int(img.size[1])
        
    final_width = min(widthsX.values())
    final_heigth = min(heigthY.values())
    print("Final Width = ")
    print(final_width)
    
    print("Final Heigth = ")
    print(final_heigth)
    
    for marker in links:
        if not os.path.exists(os.path.join(tiffdir, marker + ".tiff")):
            print(marker)
            url = os.path.join(input_dir, sample_name, links[marker])
            if os.path.isfile(url):    
                img = Image.open(url)
                 
                startx = img.size[0] - final_width
                starty = 0
                stopx = img.size[0]
                stopy = final_heigth
            
                img_cropped = img.crop((startx, starty, stopx, stopy))
                
                # Remove noise
                # background_img = img.crop(0,0,99,99)
                
                img_cropped.save(os.path.join(tiffdir, marker + '.tiff'))

    return os.path.join(output_dir, "input") 


def detect_bright_artifacts(tiff_dir, samples, markers = ["CD14", "CD16", "CD195", "HLAABC"], mask_suffix = "_spots.tiff"):
    for sample in samples:
        for marker in markers:
            url = os.path.join(tiff_dir, sample,  marker + ".tiff")
            print(marker)
            print(url)
            if os.path.isfile(url):
                image = tif.imread(url)
                image[image < 40000] = 0
                image[image > 0] = 1
                label_img = label(image)
                regions = regionprops(label_img)
                count = 0
                for i in range(len(regions)):
                    if regions[i].area < 200:
                        label_img[label_img == i+1] = 0
                    else:
                        count = count + 1
                label_img = segmentation.expand_labels(label_img, distance=70)
                print("Total spots removed for " + marker + ": " + str(count))
                tif.imwrite(os.path.join(tiff_dir, sample,  marker + mask_suffix), label_img)
            else:
                print(marker + " is not present in tiff dir.")

def remove_artifacts_from_multiple_masks(tiff_dir, samples, markers, mask_suffix = "_spots.tiff"):
    for sample in samples:
        for marker in markers:
            mask = os.path.join(tiff_dir, sample,  marker + mask_suffix)
            if os.path.isfile(mask):
                remove_artifact_from_mask(mask, tiff_dir, sample, marker)
            else:
                print("Mask for " + marker + " is not present in tiff dir.")
            
def remove_artifacts_from_single_mask(mask, tiff_dir, samples, markers):
    for sample in samples:
        for marker in markers:
            remove_artifact_from_mask(mask, tiff_dir, sample, marker)
            
def remove_artifact_from_mask(mask, tiff_dir, sample, marker):
    mask_tiff = tif.imread(mask)
    image = tif.imread(os.path.join(tiff_dir, sample,  marker + ".tiff"))
    tif.imwrite(os.path.join(tiff_dir, sample,  marker + "_original.tiff"), image)
    image[mask_tiff > 0] = 0
    tif.imwrite(os.path.join(tiff_dir, sample,  marker + ".tiff"), image)

def napari_polygon_to_mask(sample, image_shape, polygon_file, tiff_dir, mask_suffix = "_remove_mask.tiff"):
    
    poly = pd.read_csv(polygon_file)
    poly_points = np.zeros((poly.shape[0],2))
    for i in range(poly.shape[0]):
        poly_points[i] = [max(0, round(poly["axis-0"][i])), max(0, round(poly["axis-1"][i]))]
    polygon = np.array(poly_points)
    mask = polygon2mask(image_shape, polygon)
    mask = mask.astype(int)
    tif.imwrite(os.path.join(tiff_dir, sample, sample + mask_suffix), mask)

def remove_artifacts_below_threshold(thresholds, tiff_dir, sample, markers):
    for i in range(len(thresholds)):
        marker = markers[i]
        threshold = thresholds[i]
        print(marker)
        print(threshold)
        remove_artifact_below_threshold(threshold, tiff_dir, sample, marker)
            
def remove_artifact_below_threshold(threshold, tiff_dir, sample, marker):
    image = tif.imread(os.path.join(tiff_dir, sample,  marker + ".tiff"))
    tif.imwrite(os.path.join(tiff_dir, sample,  marker + "_original.tiff"), image)
    image[image < threshold] = 0
    tif.imwrite(os.path.join(tiff_dir, sample,  marker + ".tiff"), image)
    

def smooth_fill_holes(samples, markers, tiff_dir, suffix = "_filled.tiff"):
    for sample in samples:
        for marker in markers:
            image = tif.imread(os.path.join(tiff_dir, sample,  marker  + ".tiff"))
            image = gaussian(image, sigma=10)
            seed = np.copy(image)
            seed[1:-1, 1:-1] = image.max()
            mask = image
            filled = reconstruction(seed, mask, method='erosion')
            tif.imwrite(os.path.join(tiff_dir, sample, marker + suffix), mask)

###############################################################################
## General util functions
###############################################################################

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def parse_var(s):
    """
    Parse a key, value pair, separated by '='
    That's the reverse of ShellArgs.

    On the command line (argparse) a declaration will typically look like:
        foo=hello
    or
        foo="hello world"
    """
    items = s.split('=')
    key = items[0].strip() # we remove blanks around keys, as is logical
    if len(items) > 1:
        # rejoin the rest:
        value = '='.join(items[1:])
    return (key, value)



def parse_vars(items, to = "int"):
    """
    Parse a series of key-value pairs and return a dictionary
    """
    d = {}

    if items:
        for item in items:
            key, value = parse_var(item)
            if to == "int":
                d[key] = int(value)
            elif to == "float":
                d[key] = float(value)
            elif to == "bool":
                d[key] = value == 'True'
            elif to == "str":
                d[key] = str(value)
            else :
                d[key] = value
    return d


def filter_with_nuclear_mask(fovs: List, tiff_dir: str, seg_dir: str, channel: str,
                             nuc_seg_suffix: str = "_nuclear.tiff", img_sub_folder: str = None,
                             exclude: bool = True):
    """
    Filters out background staining using subcellular marker localization.

    Non-nuclear signal is removed from nuclear markers and vice-versa for membrane markers.

    Args:
        fovs (list):
            The list of fovs to filter
        tiff_dir (str):
            Name of the directory containing the tiff files
        seg_dir (str):
            Name of the directory containing the segmented files
        channel (str):
            Channel to apply filtering to
        nuc_seg_suffix (str):
            The suffix for the nuclear channel.
            (i.e. for "fov1", a suffix of "_nuclear.tiff" would make a file named
            "fov1_nuclear.tiff")
        img_sub_folder (str):
            Name of the subdirectory inside `tiff_dir` containing the tiff files.
            Set to `None` if there isn't any.
        exclude (bool):
            Whether to filter out nuclear or membrane signal
    """
    # if seg_dir is None, the user cannot run filtering
    if seg_dir is None:
        print('No seg_dir provided, you must provide one to run nuclear filtering')
        return

    # raise an error if the provided seg_dir does not exist
    io_utils.validate_paths(seg_dir)

    # convert to path-compatible format
    if img_sub_folder is None:
        img_sub_folder = ''

    for fov in fovs:
        # load the channel image in
        img = load_utils.load_imgs_from_tree(data_dir=tiff_dir, img_sub_folder=img_sub_folder,
                                             fovs=[fov], channels=[channel]).values[0, :, :, 0]

        # load the segmented image in
        seg_img_name: str = f"{fov}{nuc_seg_suffix}"
        seg_img = imread(os.path.join(seg_dir, seg_img_name))

        # mask out the nucleus
        if exclude:
            suffix = "_nuc_exclude.tiff"
            seg_mask = seg_img > 0
        # mask out the membrane
        else:
            suffix = "_nuc_include.tiff"
            seg_mask = seg_img == 0

        # filter out the nucleus or membrane depending on exclude parameter
        img[seg_mask] = 0

        # save filtered image
        image_utils.save_image(os.path.join(tiff_dir, fov, img_sub_folder, channel + suffix), img)
        return

# define the cell table path
def combine_cell_tables(cell_table_dir, samples, cell_table_prefix = "cell_table_size_normalized"):
    """
    Combine multiple single-fov matrices into a large matrix    
    """
    cell_table = []

    for samp in samples:
        df = pd.read_csv(os.path.join(cell_table_dir, samp + "_" + cell_table_prefix + ".csv.gz"))
        cell_table.append(df)
    
    cell_table = pd.concat(cell_table, axis=0)
    cell_table_path = os.path.join(cell_table_dir, cell_table_prefix + ".csv.gz")
    cell_table.to_csv(cell_table_path)
    return(cell_table_path)
    

def generate_and_save_pixel_cluster_masks(fovs: List[str],
                                          base_dir: Union[pathlib.Path, str],
                                          save_dir: Union[pathlib.Path, str],
                                          tiff_dir: Union[pathlib.Path, str],
                                          chan_file: Union[pathlib.Path, str],
                                          pixel_data_dir: Union[pathlib.Path, str],
                                          pixel_cluster_col: str = 'pixel_meta_cluster',
                                          sub_dir: str = None,
                                          name_suffix: str = ''):
    """Generates pixel cluster masks and saves them for downstream analysis.

    Args:
        fovs (List[str]):
            A list of fovs to generate and save pixel masks for.
        base_dir (Union[pathlib.Path, str]):
            The path to the data directory.
        save_dir (Union[pathlib.Path, str]):
            The directory to save the generated pixel cluster masks.
        tiff_dir (Union[pathlib.Path, str]):
            The path to the directory with the tiff data.
        chan_file (Union[pathlib.Path, str]):
            The path to the channel file inside each FOV folder (FOV folder as root).
            Used to determine dimensions of the pixel mask.
        pixel_data_dir (Union[pathlib.Path, str]):
            The path to the data with full pixel data.
            This data should also have the SOM and meta cluster labels appended.
        pixel_cluster_col (str, optional):
            The path to the data with full pixel data.
            This data should also have the SOM and meta cluster labels appended.
            Defaults to 'pixel_meta_cluster'.
        sub_dir (str, optional):
            The subdirectory to save the images in. If specified images are saved to
            `"data_dir/sub_dir"`. If `sub_dir = None` the images are saved to `"data_dir"`.
            Defaults to `None`.
        name_suffix (str, optional):
            Specify what to append at the end of every pixel mask. Defaults to `''`.
    """

    # create the pixel cluster masks across each fov
    with tqdm(total=len(fovs), desc="Pixel Cluster Mask Generation") as pixel_mask_progress:
        for fov in fovs:
            # define the path to provided channel file in the fov dir, used to calculate dimensions
            chan_file_path = os.path.join(fov, chan_file)

            # generate the pixel mask for the FOV
            pixel_mask: np.ndarray =\
                generate_pixel_cluster_mask(fov=fov, base_dir=base_dir, tiff_dir=tiff_dir,
                                            chan_file_path=chan_file_path,
                                            pixel_data_dir=pixel_data_dir,
                                            pixel_cluster_col=pixel_cluster_col)

            # save the pixel mask generated
            data_utils.save_fov_mask(fov, data_dir=save_dir, mask_data=pixel_mask, sub_dir=sub_dir,
                          name_suffix=name_suffix)

            pixel_mask_progress.update(1)


def generate_pixel_cluster_mask(fov, base_dir, tiff_dir, chan_file_path,
                                pixel_data_dir, pixel_cluster_col='pixel_meta_cluster'):
    """For a fov, create a mask labeling each pixel with their SOM or meta cluster label

    Args:
        fov (list):
            The fov to relabel
        base_dir (str):
            The path to the data directory
        tiff_dir (str):
            The path to the tiff data
        chan_file_path (str):
            The path to the sample channel file to load (`tiff_dir` as root).
            Used to determine dimensions of the pixel mask.
        pixel_data_dir (str):
            The path to the data with full pixel data.
            This data should also have the SOM and meta cluster labels appended.
        pixel_cluster_col (str):
            Whether to assign SOM or meta clusters
            needs to be `'pixel_som_cluster'` or `'pixel_meta_cluster'`

    Returns:
        numpy.ndarray:
            The image overlaid with pixel cluster labels
    """

    # path checking
    io_utils.validate_paths([tiff_dir, os.path.join(tiff_dir, chan_file_path),
                             os.path.join(base_dir, pixel_data_dir)])

    # verify the pixel_cluster_col provided is valid
    misc_utils.verify_in_list(
        provided_cluster_col=[pixel_cluster_col],
        valid_cluster_cols=['pixel_som_cluster', 'pixel_meta_cluster']
    )

    # verify the fov is valid
    misc_utils.verify_in_list(
        provided_fov_file=[fov + '.feather'],
        consensus_fov_files=os.listdir(os.path.join(base_dir, pixel_data_dir))
    )

    # read the sample channel file to determine size of pixel cluster mask
    channel_data = np.squeeze(imread(os.path.join(tiff_dir, chan_file_path)))

    # define an array to hold the overlays for the fov
    # use int16 to allow for Photoshop loading
    img_data = np.zeros((channel_data.shape[0], channel_data.shape[1]), dtype='int16')

    fov_data = feather.read_dataframe(
        os.path.join(base_dir, pixel_data_dir, fov + '.feather')
    )

    # ensure integer display and not float
    fov_data[pixel_cluster_col] = fov_data[pixel_cluster_col].astype(int)

    # get the pixel coordinates
    x_coords = fov_data['row_index'].values
    y_coords = fov_data['column_index'].values

    # convert to 1D indexing
    coordinates = x_coords * img_data.shape[1] + y_coords

    # get the cooresponding cluster labels for each pixel
    cluster_labels = list(fov_data[pixel_cluster_col])

    # assign each coordinate in pixel_cluster_mask to its respective cluster label
    img_subset = img_data.ravel()
    img_subset[coordinates] = cluster_labels
    img_data = img_subset.reshape(img_data.shape)

    return img_data