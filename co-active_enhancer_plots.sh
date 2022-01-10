#!/bin/bash 

for gene in `cat $1`
do

 Rscript /Users/chaymaeziyani/Desktop/MP_new_processing_october/task5-updated/6_enhancer_activity_plots/new_file.R $gene
 
done

