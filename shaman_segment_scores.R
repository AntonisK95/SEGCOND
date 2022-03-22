#!/usr/bin/Rscript

.libPaths("/nfs/users/tgraf/aklonizakis/R/x86_64-pc-linux-gnu-library/3.5")
library("misha")
library("shaman")


gdb.init("/users/tgraf/aklonizakis/ESCs_cells_testing_written_algorithm/Shaman_normalization/hg38")
gsetroot("/users/tgraf/aklonizakis/ESCs_cells_testing_written_algorithm/Shaman_normalization/hg38")
gdb.reload()


seg<-read.delim("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/0h/All_marks_timepoint_0h.txt_segmentation_all_chrs_0_05_breaks_1kb_windows_blaER.txt", sep="\t", header=F)
i<-Sys.getenv("var")


temporar<-seg[which(seg[,1]==i),]
final_matrix<-matrix(nrow=length(temporar[,1]),ncol=length(temporar[,1]))
final_matrix_median<-matrix(nrow=length(temporar[,1]),ncol=length(temporar[,1]))
for(k in 1:length(final_matrix[,1])){
for(j in k:length(final_matrix[,1])){

point_score<-gextract("blaer_normalized_0h", gintervals.2d(temporar[k,1], temporar[k,2], temporar[k,3], temporar[j,1], temporar[j,2], temporar[j,3]), colnames="score")[,7]

if(length(point_score>0)){
final_matrix[k,j]<-mean(point_score)
final_matrix_median[k,j]<-median(point_score)
}
else{
final_matrix[k,j]<-0
final_matrix_median[k,j]<-0

}


}
}

filename_mean<-paste("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/0h/shaman_scores_per_chromosome_per_segment/mean_normalized_hic_contacts_between_segments_of_",i,".txt",sep="")
filename_median<-paste("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/0h/shaman_scores_per_chromosome_per_segment/median_normalized_hic_contacts_between_segments_of_",i,".txt",sep="")
write.table(file=filename_mean,final_matrix,quote=F,sep="\t")
write.table(file=filename_median,final_matrix_median,quote=F,sep="\t")
