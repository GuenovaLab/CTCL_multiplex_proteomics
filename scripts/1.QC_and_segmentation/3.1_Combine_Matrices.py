#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Concatenating average expression matrices from all samples
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
import glob
argParser = argparse.ArgumentParser()
argParser.add_argument("-o", "--output-dir", help="Output directory.", type=dir_path)
args = argParser.parse_args()

print("Combining cell matrices...")

cell_table_dir = os.path.join(args.output_dir, "cell_table")

if os.path.exists(os.path.join(cell_table_dir, "cell_table_size_normalized.csv.gz")):
    os.remove(os.path.join(cell_table_dir, "cell_table_size_normalized.csv.gz"))
    

if os.path.exists(os.path.join(cell_table_dir, "cell_table_arcsinh_transformed.csv.gz")):
    os.remove(os.path.join(cell_table_dir, "cell_table_arcsinh_transformed.csv.gz"))
    
    
files =  glob.glob(os.path.join(cell_table_dir , "*cell_table_size_normalized*.csv.gz") )
samples = []
for samp in files:
    samp = os.path.basename(samp).replace("cell_table_size_normalized", "").replace("_", "").replace(".csv.gz", "")
    samples.append(samp)
print(samples)
combine_cell_tables(cell_table_dir, samples, "cell_table_size_normalized")
combine_cell_tables(cell_table_dir, samples, "cell_table_arcsinh_transformed")
