#### THIS IS THE SCRIPT USED TO ANNOTATE SEGMENTS ####
#### THE ANNOTATION IS BASED ON A ZERO-INFLATED NEGATIVE BINOMIAL MODEL ####
#### EACH SEGMENT IS ASSIGNED A P-VALUE, A LOGFC CHANGE AND AN ADJUSTED P-VALUE ####

### read the marks tables ####

marks_0h<-read.delim(sep="\t",header=T,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/All_marks_timepoint_0h.txt")
marks_1d<-read.delim(sep="\t",header=T,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/All_marks_timepoint_24h.txt")
marks_7d<-read.delim(sep="\t",header=T,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/All_marks_timepoint_7d.txt")

#### z-scale colunms ####
#### note that the z-scale function will first log2 transform the data ####

get_1z_score<-function(marks){
  backup<-marks
  backup[is.na(backup)]<-0
  marks[is.na(marks)]<-0
  means_of_log<-vector()
  sds_of_log<-vector()
  for(j in 4:ncol(marks)){
    
    means_of_log[j-3]<-mean(log2(marks[,j]+0.01))
    sds_of_log[j-3]<-sd(log2(marks[,j]+0.01))
    
  }
  
  for(i in 1:length(means_of_log)){
    
    backup[,i+3]<-(log2(backup[,i+3]+0.01)-means_of_log[i])/sds_of_log[i]
    
  }
  
  return(backup)
  
}

cooccupancy_mat_0h<-get_1z_score(marks_0h)
cooccupancy_mat_1d<-get_1z_score(marks_1d)
cooccupancy_mat_7d<-get_1z_score(marks_7d)

### transform z-scaled columns into "votes" per bin : Every bin that has a Z score of >=1 ###
### for all replicates of ATAC-seq and H3K27ac occupancy is turned into 1 while all the others at 0 ###
### changing the below function allows to manipulate which votes you think are important ###

votes_matrix<-function(cooccupancy_mat){

lines_of_interest<-which(cooccupancy_mat[,4]>=1 & cooccupancy_mat[,5]>=1 & cooccupancy_mat[,8]>=1 & cooccupancy_mat[,9]>=1) 

new_test<-cbind.data.frame(cooccupancy_mat[,c(1,2,3)],rep(0,length(cooccupancy_mat[,1])))
colnames(new_test)<-c("chr","start","end","co-occupancy")
new_test[lines_of_interest,4]<-1
return(new_test)
}
new_scores_0h<-votes_matrix(cooccupancy_mat_0h)
new_scores_1d<-votes_matrix(cooccupancy_mat_1d)
new_scores_7d<-votes_matrix(cooccupancy_mat_7d)

#### Read the new segments. Will need their length distribution + the actual coordinates to fit models and get all the metrics ###
segments_0h<-read.delim(sep="\t",header=F,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/0h/All_marks_timepoint_0h.txt_segmentation_all_chrs_0_05_breaks_1kb_windows_blaER.txt")
segments_1d<-read.delim(sep="\t",header=F,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/1d/All_marks_timepoint_24h.txt_segmentation_all_chrs_0_05_breaks_1kb_windows_blaER.txt")
segments_7d<-read.delim(sep="\t",header=F,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/0_05_segmentation_results/7d/All_marks_timepoint_7d.txt_segmentation_all_chrs_0_05_breaks_1kb_windows_blaER.txt")

#### Get a background random distribution's 
#### Need vector with segment distances and the above new_scores matrices ####

segment_length_0h<-segments_0h[,3]-segments_0h[,2]
segment_length_1d<-segments_1d[,3]-segments_1d[,2]
segment_length_7d<-segments_7d[,3]-segments_7d[,2]

randomize_segment_votes<-function(segment_length,new_scores){
mean_random_votes_per_segment_cooccupancy<-list()
for(j in 1:length(segment_length)){
  test<-segment_length[j]
  test2<-round((test)/5000)
  
  ## test with all_marks_0h and only atac-seq ##
  mean_vector<-vector()
  
  random_starts<-sample(seq(1,length(new_scores[,1])-test2,1),1000)
  random_ends<-random_starts+test2
  
  for(i in 1:length(random_starts)){
    mean_vector[i]<-sum(new_scores[random_starts[i]:random_ends[i],4])
  }
  mean_random_votes_per_segment_cooccupancy[[j]]<-mean_vector
}
return(mean_random_votes_per_segment_cooccupancy)
}

background_random_votes_0h<-randomize_segment_votes(segment_length_0h,new_scores_0h)
background_random_votes_1d<-randomize_segment_votes(segment_length_1d,new_scores_1d)
background_random_votes_7d<-randomize_segment_votes(segment_length_7d,new_scores_7d)

### Get actual vote distribution ####
get_actual_votes<-function(segments,new_scores){
actual_cooccupancy_votes<-vector()
for(i in 1:length(segments[,1])){
  
  chr<-as.character(segments[i,1])
  start<-as.numeric(segments[i,2])
  end<-as.numeric(segments[i,3])
  
  marks_start<-which(new_scores[,1]==chr & new_scores[,2]==start)
  marks_end<-which(new_scores[,1]==chr & new_scores[,3]==end)
  
  actual_cooccupancy_votes[i]<-sum(new_scores[marks_start:marks_end,4])
  
}
return(actual_cooccupancy_votes)
}
actual_votes_0h<-get_actual_votes(segments_0h,new_scores_0h)
actual_votes_1d<-get_actual_votes(segments_1d,new_scores_1d)
actual_votes_7d<-get_actual_votes(segments_7d,new_scores_7d)

#### Fit a zero-inflated negative binomial distribution to background distribution ####
#### Get a p-value and an enrichment score for each segment ####

library(fitdistrplus)
library(VGAM)
library(gamlss)  

get_zinegbin_pvalues<-function(mean_random_votes_per_segment_cooccupancy,actual_cooccupancy_votes){
size_estimates_cooccupancy<-vector()
mean_estimates_cooccupancy<-vector()
pstr0_estimates_cooccupancy<-vector()

### At certain cases the algorithm failed to find the optimal parameters : Tweaking the starting parameters seemed to resolve the problem, hence the multiple if statements below ###

for(i in 1:length(mean_random_votes_per_segment_cooccupancy)){
  fit<-"NA"
  check<-try(fit <- fitdist(mean_random_votes_per_segment_cooccupancy[[i]],lower = c(0,0,0), "zinegbin",start = list(pstr0=(length(which(mean_random_votes_per_segment_cooccupancy[[i]]==0))/length(mean_random_votes_per_segment_cooccupancy[[i]])),size=(mean(mean_random_votes_per_segment_cooccupancy[[i]])^2)/(var(mean_random_votes_per_segment_cooccupancy[[i]])-mean(mean_random_votes_per_segment_cooccupancy[[i]])),munb = mean(mean_random_votes_per_segment_cooccupancy[[i]]))))
  
  if(class(check)=="try-error"){
    check<-try(fit <- fitdist(mean_random_votes_per_segment_cooccupancy[[i]], "zinegbin",start = list(pstr0=0.8,size=(mean(mean_random_votes_per_segment_cooccupancy[[i]])^2)/(var(mean_random_votes_per_segment_cooccupancy[[i]])-mean(mean_random_votes_per_segment_cooccupancy[[i]])),munb = mean(mean_random_votes_per_segment_cooccupancy[[i]]))))
    
    if(class(check)=="try-error"){
      
      fit <- try(fitdist(mean_random_votes_per_segment_cooccupancy[[i]], "zinegbin",start = list(pstr0=0.99,size=(mean(mean_random_votes_per_segment_cooccupancy[[i]])^2)/(var(mean_random_votes_per_segment_cooccupancy[[i]])-mean(mean_random_votes_per_segment_cooccupancy[[i]])),munb = mean(mean_random_votes_per_segment_cooccupancy[[i]]))))
      
      if(class(check)=="try-error"){
        fit <- fitdist(mean_random_votes_per_segment_cooccupancy[[i]], "zinegbin",start = list(pstr0=0.8,size=0.6,munb = 0.6))
        pstr0_estimates_cooccupancy[i]<-as.numeric(fit$estimate[1])
        size_estimates_cooccupancy[i]<-as.numeric(fit$estimate[2])
        mean_estimates_cooccupancy[i]<-as.numeric(fit$estimate[3])
        
      }
      
      else{
      pstr0_estimates_cooccupancy[i]<-as.numeric(fit$estimate[1])
      size_estimates_cooccupancy[i]<-as.numeric(fit$estimate[2])
      mean_estimates_cooccupancy[i]<-as.numeric(fit$estimate[3])
    }}
    else{
    pstr0_estimates_cooccupancy[i]<-as.numeric(fit$estimate[1])
    size_estimates_cooccupancy[i]<-as.numeric(fit$estimate[2])
    mean_estimates_cooccupancy[i]<-as.numeric(fit$estimate[3])
    }
  }
  else{
    pstr0_estimates_cooccupancy[i]<-as.numeric(fit$estimate[1])
    size_estimates_cooccupancy[i]<-as.numeric(fit$estimate[2])
    mean_estimates_cooccupancy[i]<-as.numeric(fit$estimate[3])
  }
  
}
pvalues_cooccupancy<-vector();enrichment_cooccupancy<-vector()
for(i in 1:length(actual_cooccupancy_votes)){
  
  pvalues_cooccupancy[i]<-1-pzinegbin(actual_cooccupancy_votes[i],size=size_estimates_cooccupancy[i],munb=mean_estimates_cooccupancy[i],pstr0=pstr0_estimates_cooccupancy[i])
  enrichment_cooccupancy[i]<-log2((actual_cooccupancy_votes[i]+0.01)/(mean_estimates_cooccupancy[i]+0.01))
  
}

return(list(enrichment_cooccupancy,pvalues_cooccupancy))
}

metrics_0h<-get_zinegbin_pvalues(background_random_votes_0h,actual_votes_0h)
metrics_1d<-get_zinegbin_pvalues(background_random_votes_1d,actual_votes_1d)
metrics_7d<-get_zinegbin_pvalues(background_random_votes_7d,actual_votes_7d)


### attach a table with segments and logfc changes + p-values ###

final_0h<-cbind.data.frame(segments_0h,metrics_0h[[1]],metrics_0h[[2]]);colnames(final_0h)<-c("chr","start","end","log2FC","p-value")
final_1d<-cbind.data.frame(segments_1d,metrics_1d[[1]],metrics_1d[[2]]);colnames(final_1d)<-c("chr","start","end","log2FC","p-value")
final_7d<-cbind.data.frame(segments_7d,metrics_7d[[1]],metrics_7d[[2]]);colnames(final_7d)<-c("chr","start","end","log2FC","p-value")

### write tables and terminate script ###

write.table(final_0h,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/New_annotation_files/Segments_0d_with_zinegbin_metrics.txt",quote=F,row.names=F,sep="\t")
write.table(final_1d,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/New_annotation_files/Segments_1d_with_zinegbin_metrics.txt",quote=F,row.names=F,sep="\t")
write.table(final_7d,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/New_annotation_files/Segments_7d_with_zinegbin_metrics.txt",quote=F,row.names=F,sep="\t")
 

##### get pvalues/logfcs for segments of timepoint X with respect to all other timepoints! ######
#### 3x3 operations ####


background_random_votes_0h_1d<-randomize_segment_votes(segment_length_0h,new_scores_1d)
background_random_votes_0h_7d<-randomize_segment_votes(segment_length_0h,new_scores_7d)

background_random_votes_1d_0h<-randomize_segment_votes(segment_length_1d,new_scores_0h)
background_random_votes_1d_7d<-randomize_segment_votes(segment_length_1d,new_scores_7d)


background_random_votes_7d_0h<-randomize_segment_votes(segment_length_7d,new_scores_0h)
background_random_votes_7d_1d<-randomize_segment_votes(segment_length_7d,new_scores_1d)


### Get actual vote distribution ####
get_actual_votes<-function(segments,new_scores){
  actual_cooccupancy_votes<-vector()
  for(i in 1:length(segments[,1])){
    
    chr<-as.character(segments[i,1])
    start<-as.numeric(segments[i,2])
    end<-as.numeric(segments[i,3])
    
    marks_start<-which(new_scores[,1]==chr & new_scores[,2]==start)
    marks_end<-which(new_scores[,1]==chr & new_scores[,3]==end)
    
    actual_cooccupancy_votes[i]<-sum(new_scores[marks_start:marks_end,4])
    
  }
  return(actual_cooccupancy_votes)
}

actual_votes_0h_1d<-get_actual_votes(segments_0h,new_scores_1d)
actual_votes_0h_7d<-get_actual_votes(segments_0h,new_scores_7d)

actual_votes_1d_0h<-get_actual_votes(segments_1d,new_scores_0h)
actual_votes_1d_7d<-get_actual_votes(segments_1d,new_scores_7d)

actual_votes_7d_0h<-get_actual_votes(segments_7d,new_scores_0h)
actual_votes_7d_1d<-get_actual_votes(segments_7d,new_scores_1d)


metrics_0h_1d<-get_zinegbin_pvalues(background_random_votes_0h_1d,actual_votes_0h_1d)
metrics_0h_7d<-get_zinegbin_pvalues(background_random_votes_0h_7d,actual_votes_0h_7d)
metrics_1d_0h<-get_zinegbin_pvalues(background_random_votes_1d_0h,actual_votes_1d_0h)
metrics_1d_7d<-get_zinegbin_pvalues(background_random_votes_1d_7d,actual_votes_1d_7d)

metrics_7d_0h<-get_zinegbin_pvalues(background_random_votes_7d_0h,actual_votes_7d_0h)
metrics_7d_1d<-get_zinegbin_pvalues(background_random_votes_7d_1d,actual_votes_7d_1d)

### merge for each timepoint all the info together ###


final_0h_all<-cbind.data.frame(segments_0h,metrics_0h[[1]],metrics_0h[[2]],metrics_0h_1d[[1]],metrics_0h_1d[[2]],metrics_0h_7d[[1]],metrics_0h_7d[[2]]);colnames(final_0h_all)<-c("chr","start","end","log2FC_0h","p-value_0h","log2FC_1d","p-value_1d","log2FC_7d","p-value_7d")
final_1d_all<-cbind.data.frame(segments_1d,metrics_1d_0h[[1]],metrics_1d_0h[[2]],metrics_1d[[1]],metrics_1d[[2]],metrics_1d_7d[[1]],metrics_1d_7d[[2]]);colnames(final_1d_all)<-c("chr","start","end","log2FC_0h","p-value_0h","log2FC_1d","p-value_1d","log2FC_7d","p-value_7d")
final_7d_all<-cbind.data.frame(segments_7d,metrics_7d_0h[[1]],metrics_7d_0h[[2]],metrics_7d_1d[[1]],metrics_7d_1d[[2]],metrics_7d[[1]],metrics_7d[[2]]);colnames(final_7d_all)<-c("chr","start","end","log2FC_0h","p-value_0h","log2FC_1d","p-value_1d","log2FC_7d","p-value_7d")


write.table(final_0h_all,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/New_annotation_files/Segments_0d_with_zinegbin_metrics_for_all_timepoints.txt",quote=F,row.names=F,sep="\t")
write.table(final_1d_all,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/New_annotation_files/Segments_1d_with_zinegbin_metrics_for_all_timepoints.txt",quote=F,row.names=F,sep="\t")
write.table(final_7d_all,file="/users/tgraf/aklonizakis/blaer_cells_segmentation/New_annotation_files/Segments_7d_with_zinegbin_metrics_for_all_timepoints.txt",quote=F,row.names=F,sep="\t")





