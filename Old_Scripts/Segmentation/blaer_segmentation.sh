#!/usr/bin/bash
#$ -q long-sl7
#$ -l h_rt=24:0:0
#$ -l virtual_free=10G

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load R/3.5.0-foss-2018a-X11-20180131

myarr=(/users/tgraf/aklonizakis/blaer_cells_segmentation/All_marks_timepoint_0h.txt /users/tgraf/aklonizakis/blaer_cells_segmentation/All_marks_timepoint_24h.txt /users/tgraf/aklonizakis/blaer_cells_segmentation/All_marks_timepoint_7d.txt)

export var="${myarr[${SGE_TASK_ID}-1]}"

Rscript /users/tgraf/aklonizakis/blaer_cells_segmentation/blaer_segmentation.R


