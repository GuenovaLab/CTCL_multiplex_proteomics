#!/bin/bash

# Running celltype classication model evaluation using CellSighter
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

samples=$(ls /mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/output/input | grep ROI | tr "\n" " ")

samples="ROI-06-Annotator1 ROI-07-Annotator1 ROI-15-Annotator3 ROI-15-Annotator1 ROI-17-Annotator3" 


rootdir=/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/output/CellSighter/celltype/Predictions/

for sample in $(echo $samples);
	
	    do




		echo "###############################################################################"
		echo "###############################################################################"
        	echo $sample
		echo "###############################################################################"
		echo "###############################################################################"
	
		
		cp $rootdir/config_template.json $rootdir/config.json


		sed -i "s+AAAA+\"$sample\"+g" $rootdir/config.json


	       	ulimit -n 100000000; python /mnt/RECHERCHE/GUENOVA_LAB/Pacome/GitLab/CellSighter-main/eval.py --base_path=$rootdir
		mv $rootdir/val_results.csv $rootdir/${sample}_val_results.csv

		echo
		echo
		echo



done;


