#!/bin/bash


# Running positive vs negative marker classication model evaluation using CellSighter
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

markers=$(cat /mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/annotation/channels.txt | tr "\n" " ")
samples=$(ls /mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/output/input | grep ROI | tr "\n" " ")


for marker in $(echo markers);
	
	do
		trainingdir=/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/output/CellSighter/marker/marker_classification/$marker
		markerdir=/mnt/RECHERCHE/GUENOVA_LAB/Project_Multiplex_Phenotyping/Miltenyi/Muliplex_Imaging_Pipeline/output/CellSighter/marker/Predictions/$marker
		channels=$markerdir/channels.txt
	
		
		for validation in $(ls $trainingdir/val*.csv | tr "\n" " ");
		do

		   iteration=$(basename $validation | sed 's/val_results_\|.csv//g')

		   
		   awk -v iter=$iteration -v FS="," 'NR>1{if($2==0 && $4==0){TN+=1}; if($2==0 && $4==1){FN+=1}; if($2==1 && $4==0){FP+=1}; if($2==1 && $4==1){TP+=1};} \
		   END{print iter" "TP / (TP + 0.5 *  (FP + FN))}' $validation	
		   
		done | sort -n -k2 | tail -1 | cut -f1 -d" " > $trainingdir/best_iter.txt
	
	
		weights=$trainingdir/weights_$(cat $trainingdir/best_iter.txt)_count.pth
	    
		head -n1 $markerdir/../val_results.csv > $markerdir/${marker}_predictions.csv
		    
		    for sample in $(echo $samples);
		    do
			echo "###############################################################################"
			echo "###############################################################################"
 	        	echo $marker - $sample
			echo "###############################################################################"
			echo "###############################################################################"
		
			
			if [[ -f $markerdir/CellTypes/cells/${sample}.npz ]]
			then


		       	ulimit -n 1000000; python /mnt/RECHERCHE/GUENOVA_LAB/Pacome/GitLab/CellSighter-main/eval.py --base_path=$markerdir
			mv $markerdir/val_results.csv $markerdir/${sample}_val_results.csv
			tail -n +2 $markerdir/${sample}_val_results.csv >> $markerdir/${marker}_predictions.csv
			echo
			echo
			echo
			fi
	done;


done;


