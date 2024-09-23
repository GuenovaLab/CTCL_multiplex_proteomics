#!/bin/bash

# Running quality control, artifact removal, segmentation and 
# pixel classification
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

input_dir=$1
output_dir=$2

mkdir -p $output_dir
samples=$(ls $input_dir | grep ROI | sed 's/\///g' | tr "\n" ",")

export PATH="$PATH":"/usr/local/cuda-12.2/bin/"
export LD_LIBRARY_PATH="/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cublas/lib/":\
"/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cufft/lib/":\
"/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/curand/lib/":\
"/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cusolver/lib/":\
"/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cusparse/lib/":\
"/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cudnn/lib/":\
"/home/localadmin/anaconda3/lib/python3.10/site-packages/nvidia/cuda_runtime/lib/":"$LD_LIBRARY_PATH"

./scripts/1.QC_and_segmentation/2.0_Input_QC_and_filtering.py -i $input_dir/ -o $output_dir -s $samples -m annotation/marker_metadata.csv -a annotation/

samples=$(ls $input_dir | grep ROI | sed 's/\///g' | tr "\n" " ")

for sample_name in  $(echo $samples  | tr "," " ")
do
	echo
	echo $sample_name
        echo ./scripts/1.QC_and_segmentation/3.0_DeepCell_Nuclear_Detection.py -i $input_dir/ -o $output_dir -s ${sample_name} -n DAPI -c Cytokeratin,CD45RO,betaActin,Galectin3,p16INK4a -m annotation/marker_metadata.csv --resolution 0.5 -C interior_threshold=0.1 maxima_threshold=0.5 -N interior_threshold=0.05 maxima_threshold=0.5
done
	
./scripts/1.QC_and_segmentation/3.1_Combine_Matrices.py -o $output_dir

# Run using jupyter notebook "4.0_Pixel_Classification.ipynb"
# ./scripts/1.QC_and_segmentation/4.0_Pixel_Classification.py -i $input_dir -o $output_dir -n All -p 0.05 -f CD3="True" CD4="True" CD8a="True" CD45RA="True" CD45RO="True" HLADR="True" FoxP3="True" CD5="True" CD196="True" CD204="True" HLADRDPDQ="True" CD13="True" Cytokeratin="True" CD14="True"  -m "CD3,CD4,CD8a,FoxP3" 

./scripts/1.QC_and_segmentation/4.1_Pixel_Cell_Classification.py -o $output_dir -n region -k 12




