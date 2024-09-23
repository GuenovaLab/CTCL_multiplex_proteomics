#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pixel Cell Classification using Pixie
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

from MultiplexImaging_utils import *



argParser = argparse.ArgumentParser()
argParser.add_argument("-o", "--output-dir", help="Output Directory", type=dir_path)
argParser.add_argument("-n", "--name", help="Prefix to add for the analysis.", type=str, default = "example")
argParser.add_argument("-k", "--max-k", help="The number of consensus clusters desired. Defaults to 20.", type=int, default=20)
args = argParser.parse_args()


###############################################################################
## Loading
###############################################################################

print()
print("Running Cell Pixel Classification...")
print()
print("######################################################################################")     
print("Output directory = %s" % args.output_dir)
print("Name = %s" % args.name)
samples = os.listdir(os.path.join(args.output_dir, "input"))

print("Samples = %s" % samples)
print("Maximum K = %s" % args.max_k)

segmentation_dir = os.path.join(args.output_dir, "segmentation")
print("Seg Dir " )
print(segmentation_dir)
cell_table_dir = os.path.join(args.output_dir, "cell_table")
tiff_dir = os.path.join(args.output_dir, "input")
# define the name of the folder containing the pixel cluster data
pixel_cluster_prefix = args.name
pixel_output_dir = os.path.join("pixie", "%s_pixel_output_dir" % pixel_cluster_prefix)
# define the name of the cell clustering params file
cell_clustering_params_name = 'cell_clustering_params.json'

print("######################################################################################")    
print()

# import required packages
import json
import os
import subprocess
from datetime import datetime as dt

import feather
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import seaborn as sns
import xarray as xr
from matplotlib import rc_file_defaults
from alpineer import io_utils, load_utils

from ark.analysis import visualize
from ark.phenotyping import cell_cluster_utils
from ark.utils import data_utils, example_dataset, plot_utils
from ark.utils.metacluster_remap_gui import (MetaClusterData, MetaClusterGui,
                                             colormap_helper,
                                             metaclusterdata_from_files)
import time
from ark.phenotyping.weighted_channel_comp import compute_p2c_weighted_channel_avg
from ark.phenotyping.cell_som_clustering import train_cell_som, cluster_cells, generate_som_avg_files 
from ark.phenotyping.cell_meta_clustering import cell_consensus_cluster, generate_meta_avg_files
from ark.phenotyping.weighted_channel_comp import generate_wc_avg_files
from ark.phenotyping.weighted_channel_comp import generate_weighted_channel_avg_heatmap

start_time = time.time()
# define the name of the cell clustering params file
cell_clustering_params_name = 'cell_clustering_params.json'

# load the params
with open(os.path.join(args.output_dir,  pixel_output_dir, cell_clustering_params_name)) as fh:
    cell_clustering_params = json.load(fh)
    
# assign the params to variables
fovs = cell_clustering_params['fovs']
channels = cell_clustering_params['channels']
segmentation_dir = cell_clustering_params['segmentation_dir']
seg_suffix = cell_clustering_params['seg_suffix']
pixel_data_dir = cell_clustering_params['pixel_data_dir']
pc_chan_avg_som_cluster_name = cell_clustering_params['pc_chan_avg_som_cluster_name']
pc_chan_avg_meta_cluster_name = cell_clustering_params['pc_chan_avg_meta_cluster_name']

# define the cell table path
cell_table_path = os.path.join(cell_table_dir, 'cell_table_size_normalized.csv.gz')



# explicitly set cell_cluster_prefix to override datetime default
cell_cluster_prefix = args.name

# define the base output cell folder
cell_output_dir = '%s_cell_output_dir' % cell_cluster_prefix
if not os.path.exists(os.path.join(args.output_dir, "pixie", cell_output_dir)):
    os.mkdir(os.path.join(args.output_dir, "pixie", cell_output_dir))
    
# define the paths to cell clustering files, explicitly set the variables to use custom names
cell_som_weights_name = os.path.join("pixie", cell_output_dir, 'cell_som_weights.feather')
cluster_counts_name = os.path.join("pixie", cell_output_dir, 'cluster_counts.feather')
cluster_counts_size_norm_name = os.path.join("pixie", cell_output_dir, 'cluster_counts_size_norm.feather')
weighted_cell_channel_name = os.path.join("pixie", cell_output_dir, 'weighted_cell_channel.feather')
cell_som_cluster_count_avg_name = os.path.join("pixie", cell_output_dir, 'cell_som_cluster_count_avg.csv')
cell_meta_cluster_count_avg_name = os.path.join("pixie", cell_output_dir, 'cell_meta_cluster_count_avg.csv')
cell_som_cluster_channel_avg_name = os.path.join("pixie", cell_output_dir, 'cell_som_cluster_channel_avg.csv')
cell_meta_cluster_channel_avg_name = os.path.join("pixie", cell_output_dir, 'cell_meta_cluster_channel_avg.csv')
cell_meta_cluster_remap_name = os.path.join("pixie", cell_output_dir, 'cell_meta_cluster_mapping.csv')    

###############################################################################
## Assigning Pixels to Clusters
###############################################################################

print("Cell to pixel clusters")
    
# define the type of pixel cluster to aggregate on
pixel_cluster_col = 'pixel_meta_cluster_rename'

# depending on which pixel_cluster_col is selected, choose the pixel channel average table accordingly
if pixel_cluster_col == 'pixel_som_cluster':
    pc_chan_avg_name = pc_chan_avg_som_cluster_name
elif pixel_cluster_col == 'pixel_meta_cluster_rename':
    pc_chan_avg_name = pc_chan_avg_meta_cluster_name
    
# generate the preprocessed data before 
cluster_counts, cluster_counts_size_norm = cell_cluster_utils.create_c2pc_data(
    samples, os.path.join(args.output_dir, pixel_data_dir), cell_table_path, pixel_cluster_col
)

# define the count columns found in cluster_counts_norm
cell_som_cluster_cols = cluster_counts_size_norm.filter(
    regex=f'{pixel_cluster_col}.*'
).columns.values

# write the unnormalized input data to cluster_counts_name for reference
feather.write_dataframe(
    cluster_counts,
    os.path.join(args.output_dir, cluster_counts_name),
    compression='uncompressed'
)



# generate the weighted cell channel expression data
pixel_channel_avg = pd.read_csv(os.path.join(args.output_dir, pc_chan_avg_name))
pixel_channel_avg["pixel_meta_cluster_rename"] = pixel_channel_avg["pixel_meta_cluster"]
weighted_cell_channel = compute_p2c_weighted_channel_avg(
    pixel_channel_avg,
    channels,
    cluster_counts,
    fovs=samples,
    pixel_cluster_col= "pixel_meta_cluster_rename"
)

# write the data to weighted_cell_channel_name
feather.write_dataframe(
    weighted_cell_channel,
    os.path.join(args.output_dir, weighted_cell_channel_name),
    compression='uncompressed'
)
    
###############################################################################
## Cell-level Pixel Clustering
###############################################################################

print("Training Cell-Level SOM Pixel clustering")
# create the cell-level SOM weights
cell_pysom = train_cell_som(
    samples,
    args.output_dir,
    cell_table_path=cell_table_path,
    cell_som_cluster_cols=cell_som_cluster_cols,
    cell_som_input_data=cluster_counts_size_norm,
    som_weights_name=cell_som_weights_name,
    num_passes=1,
    seed=42
)

# use cell SOM weights to assign cell clusters
cluster_counts_size_norm = cluster_cells(
    args.output_dir,
    cell_pysom,
    cell_som_cluster_cols=cell_som_cluster_cols
)

# generate the SOM cluster summary files
generate_som_avg_files(
    args.output_dir,
    cluster_counts_size_norm,
    cell_som_cluster_cols=cell_som_cluster_cols,
    cell_som_expr_col_avg_name=cell_som_cluster_count_avg_name
)


print("Run pixel consensus clustering")
cap = 3

# run hierarchical clustering based on cell SOM cluster assignments
cell_cc, cluster_counts_size_norm = cell_consensus_cluster(
    args.output_dir,
    cell_som_cluster_cols=cell_som_cluster_cols,
    cell_som_input_data=cluster_counts_size_norm,
    cell_som_expr_col_avg_name=cell_som_cluster_count_avg_name,
    max_k=args.max_k,
    cap=cap
)

# generate the meta cluster summary files
generate_meta_avg_files(
    args.output_dir,
    cell_cc,
    cell_som_cluster_cols=cell_som_cluster_cols,
    cell_som_input_data=cluster_counts_size_norm,
    cell_som_expr_col_avg_name=cell_som_cluster_count_avg_name,
    cell_meta_expr_col_avg_name=cell_meta_cluster_count_avg_name
)

# generate weighted channel summary files
generate_wc_avg_files(
    samples,
    channels,
    args.output_dir,
    cell_cc,
    cell_som_input_data=cluster_counts_size_norm,
    weighted_cell_channel_name=weighted_cell_channel_name,
    cell_som_cluster_channel_avg_name=cell_som_cluster_channel_avg_name,
    cell_meta_cluster_channel_avg_name=cell_meta_cluster_channel_avg_name
)
    
###############################################################################
## Exporting
###############################################################################

print("Plotting Heatmap")
rc_file_defaults()
plt.ion()

cell_mcd = metaclusterdata_from_files(
    os.path.join(args.output_dir, cell_som_cluster_count_avg_name),
    cluster_type='cell',
    prefix_trim=pixel_cluster_col + '_'
)
cell_mcd.output_mapping_filename = os.path.join(args.output_dir, cell_meta_cluster_remap_name)
cell_mcg = MetaClusterGui(cell_mcd, width=10)


# df = pd.DataFrame(cell_mcd.mapping)
# df["mc_name"] = df['metacluster']
# df.rename(columns={"cluster": "pixel_som_cluster", "metacluster": "pixel_meta_cluster", "mc_name": "pixel_meta_cluster_rename"})
# df.to_csv( os.path.join(args.output_dir, cell_meta_cluster_remap_name))

raw_cmap, renamed_cmap = colormap_helper.generate_meta_cluster_colormap_dict(
    meta_cluster_remap_path = os.path.join(args.output_dir, cell_meta_cluster_remap_name),
    cmap = cell_mcg.im_cl.cmap
)

df = pd.read_csv(os.path.join(args.output_dir, cell_som_cluster_channel_avg_name))
df["cell_meta_cluster_rename"] = df["cell_meta_cluster"]
df.to_csv(os.path.join(args.output_dir, cell_som_cluster_channel_avg_name))

generate_weighted_channel_avg_heatmap(
    os.path.join(args.output_dir, cell_som_cluster_channel_avg_name),
    'cell_som_cluster',
    channels,
    raw_cmap,
    renamed_cmap
)
plt.savefig(os.path.join(args.output_dir, "pixie", cell_output_dir, 'heatmap_weighted_channel_avg.pdf'))

df = pd.read_csv(os.path.join(args.output_dir, cell_meta_cluster_channel_avg_name))
df["cell_meta_cluster_rename"] = df["cell_meta_cluster"]
df.to_csv(os.path.join(args.output_dir, cell_meta_cluster_channel_avg_name))

generate_weighted_channel_avg_heatmap(
    os.path.join(args.output_dir, cell_meta_cluster_channel_avg_name),
    'cell_meta_cluster_rename',
    channels,
    raw_cmap,
    renamed_cmap
)
plt.savefig(os.path.join(args.output_dir, "pixie", cell_output_dir, 'heatmap_weighted_channel_avg_metaclusters.pdf'))

###############################################################################
## Exporting Restults
###############################################################################
subset_cell_fovs = samples

# generate and save the cell  cluster masks for each fov in subset_cell_fovs
data_utils.generate_and_save_cell_cluster_masks(
    fovs=subset_cell_fovs,
    base_dir=args.output_dir,
    save_dir=os.path.join(args.output_dir, "pixie", cell_output_dir),
    seg_dir=segmentation_dir,
    cell_data=cluster_counts_size_norm,
    seg_suffix=seg_suffix,
    sub_dir='cell_masks',
    name_suffix='_cell_mask'
)


for cell_fov in subset_cell_fovs:
    cell_cluster_mask = load_utils.load_imgs_from_dir(
        data_dir = os.path.join(args.output_dir, "pixie", cell_output_dir, "cell_masks"),
        files=[cell_fov + "_cell_mask.tiff"],
        trim_suffix="_cell_mask",
        match_substring="_cell_mask",
        xr_dim_name="cell_mask",
        xr_channel_names=None,
    )

    plot_utils.plot_pixel_cell_cluster_overlay(
        cell_cluster_mask,
        [cell_fov],
        os.path.join(args.output_dir, cell_meta_cluster_remap_name),
        metacluster_colors=raw_cmap
    )
    
    plt.savefig(os.path.join(args.output_dir, "pixie", cell_output_dir, "cell_masks", cell_fov + "_overlay_meta_clusters" +'.pdf'))

cluster_counts_size_norm["cell_meta_cluster_rename"] = cluster_counts_size_norm['cell_meta_cluster']
cell_cluster_utils.add_consensus_labels_cell_table(
    args.output_dir, cell_table_path, cluster_counts_size_norm
)

feather.write_dataframe(
    cluster_counts_size_norm,
    os.path.join(args.output_dir, cluster_counts_size_norm_name),
    compression='uncompressed'
)
print('Finished Running Cell Pixel Classification in %s seconds', round(time.time() - start_time, 2))
