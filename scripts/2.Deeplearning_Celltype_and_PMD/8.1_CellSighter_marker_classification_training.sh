#!/bin/bash

# Running positive vs negative marker classication model training using CellSighter
# authors: Pacome Prompsy
# contact: pacome.prompsy@chuv.ch
# Guenova Lab
# CHUV (Centre Hospitalier Universitaire Vaudois), Lausanne, Suisse

# CellSighter:
# https://github.com/KerenLab/CellSighter
#
# Amitay, Y., Bussi, Y., Feinstein, B. et al. CellSighter: a neural network 
# to classify cells in highly multiplexed images. 
# Nat Commun 14, 4302 (2023). https://doi.org/10.1038/s41467-023-40066-7

markers=$(cat annotation/channels_same_range.txt | tr "\n" " ")

for marker in $(echo markers)
	    do
	   
		echo "###############################################################################"
		echo "###############################################################################"
 	        echo $marker
		echo "###############################################################################"
		echo "###############################################################################"

		python ../../../Pacome/GitLab/CellSighter-main/train.py --base_path=output/CellSighter/marker/marker_classification/$marker
	        
		echo
		echo
		echo
	    done;


