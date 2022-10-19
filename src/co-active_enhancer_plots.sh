#!/bin/bash 

for gene in `cat $1`
do

 Rscript Rscript.R $gene
 
done

