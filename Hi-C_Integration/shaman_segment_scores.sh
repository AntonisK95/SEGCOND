#!/usr/bin/bash

#$ -l virtual_free=40G
#$ -q long-sl7
#$ -l h_rt=72:0:0


module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load R/3.5.0-foss-2018a-X11-20180131
### all chromosomes ####

myarr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX")

export var="${myarr[${SGE_TASK_ID}-1]}"

Rscript /users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/0h/shaman_segment_scores.R
