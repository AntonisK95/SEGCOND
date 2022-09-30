#!/usr/bin/Rscript

library(strucchange)
file<-Sys.getenv("var")
temporar<-strsplit(file,split="/");name_string<-temporar[[1]][length(temporar[[1]])]
#### LOAD DATA . Format (tab separated): chr start end value1 ... value n #####

data<-read.delim(file, header=T, sep="\t")
data[is.na(data)]<-0
### REDUCE DIMENSIONS OF DATASETS TO A SINGLE VARIABLE ###

test<-prcomp(data[,4:ncol(data)],scale=T)

### IMPORTANT STEP : CHECK THE AMOUNT OF VARIANCE CAPTURED BY THE FIRST PRINCIPLE COMPONENT ###

pca_res<-summary(test)

### PREPARE DATA ###

final<-test$x[,1]
final<-cbind.data.frame(data[,1:3],final)
data<-final
colnames(data)<-c("chrom", "start", "end", "pca1")

#### NEEDED FUNCTIONS FOR THE RUNNING OF THE ALGORITHM ###

### MERGING SIMILAR BREAKS : IF A BREAK IS 5 POSITIONS NEXT TO ANOTHER, THE SMALLEST IS TO BE KEPT ###
clearence<-function(input_breaks){
  j<-1
  
  backup<-input_breaks
  while(j<=length(backup)-1){
    
    check<-backup[j+1]-backup[j]
    if(check>=5){
      j<-j+1  
      
    }
    else{
      
      backup<-backup[-(j+1)]
    }
  }
  return(backup)
}

#### DRAW FINAL COORDINATES OF SEGMENTS BASED ON THE GENERATED BREAKS ####

merge_segments<-function(choosen){
  chrom<-ch
  segments<-data.frame()
  for(i in 0:length(choosen)){
    
    if(i==0){
      
      temp1<-data1[1,2]
      temp2<-data1[choosen[1],3]
      variable<-cbind.data.frame(chrom,temp1,temp2)
      segments<-rbind.data.frame(segments,variable)
      
      
      
    }
    
    else if(i==length(choosen)){
      temp1<-data1[choosen[i]+1,2]
      temp2<-data1[size,3]
      variable<-cbind.data.frame(chrom,temp1,temp2)
      segments<-rbind.data.frame(segments,variable)
      
    }
    else{
      
      temp1<-data1[choosen[i]+1,2]
      temp2<-data1[choosen[i+1],3]
      variable<-cbind.data.frame(chrom,temp1,temp2)
      segments<-rbind.data.frame(segments,variable)
      
      
    }
  }
  return(segments)
}

##### THIS FUNCTION CALCLUATES THE MEAN AND SD OF PCA VALUES WITHIN EACH SEGMENT ###

variability_domains<-function(domains){
  sd_vector<-c()
  mean_vector<-c()
  for(i in 1:length(domains[,1])){
    line_start<-which(data1[,2]==domains[i,2])
    line_end<-which(data1[,3]==domains[i,3])
    
    values<-data1[line_start:line_end,4]
    
    sd_vector[i]<-sd(values)
    mean_vector[i]<-mean(values)
    
    
  }
  return(cbind.data.frame(mean_vector,sd_vector))
}

### LISTS TO BE FILLED WITH RESULTS ###
l<-1
segment_list<-list();testing_list<-list();breaks_list<-list()

### PREPARE CHROMOSOMES TO BE TESTED ###
### NOTE: THIS IS A MOUSE CHROMOSOME VECTOR ###

chr_list<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

### INITIATE MAIN SCRIPT ###

for(ch in chr_list){
  grep(paste(ch,"$",sep=""), data[,1])->i;data1<-data[i,c(1,2,3,4)]
  
  
  ### CREATE TIME SERIES OBJECT ###
  tsfc<-ts(data1[,4])
  size<-length(data1[,1])
  
  ##### DECIDE THE SIZE OF A SLIDING WINDOW : TESTED 500, 750, 1000 BINS AND 1000 BINS SEEMED TO OPERATE WELL ####
  
  ends<-seq(1000,size,500)
  starts<-seq(1,size-1000,500)
  final_box<-c(size-1000,size)
  breaks<-c()
  
  ### IMPORTANT : THE h PARAMETER IS THE MINIMUM SIZE A SEGMENT SHOULD HAVE. THIS IS A NON TRIVIAL PARAMETER ###
  ### THAT HAS TO BE EVALUATED SOMEHOW. AFTER TESTING, A 5% VALUE SEEMED APPROPRIATE ####
  
  for(i in 1:length(ends)){
    testing<-tsfc[starts[i]:ends[i]]
    breakpoints(testing ~ 1, h=0.05)->bp;
    results<-bp$breakpoints
    results<-results+(starts[i]-1)
    breaks<-c(breaks,results)
    
  }
  testing<-tsfc[final_box[1]:final_box[2]]
  breakpoints(testing ~ 1, h=0.05)->bp;
  results<-bp$breakpoints
  results<-results+(final_box[1]-1)
  breaks<-c(breaks,results)
  
  #### THE BREAKS VECTOR CONTAINS NON FILTERED BORDERS CALCULATED BY THE SLIDING WINDOW APPROACH ####
  
  breaks<-sort(unique(breaks))
  
  #### MERGE BREAKS THAT ARE VERY CLOSE TO EACH OTHER / IDENTICAL 
  breaks<-clearence(breaks)
  
  ### VISUALIZNG DIFFERENT BREAKS ###
  title<-paste("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/",name_string,"_",ch,"_0_05_segmentation.pdf",sep="")
  pdf(file=title, width=10)
  plot(tsfc,main=ch,ylab="PC1 smoothed values",xlab="Chromosome position")
  abline(lwd=2,col="steelblue 4",v=breaks)
  dev.off()
  
  #### FILLING LISTS WITH RESULTS ! ####
  
  final<-merge_segments(breaks)
  testing<-variability_domains(final)
  segment_list[[l]]<-final
  testing_list[[l]]<-testing
  breaks_list[[l]]<-breaks
  l<-l+1
  
}


names(segment_list)<-chr_list
names(testing_list)<-chr_list
names(breaks_list)<-chr_list

final_list<-list(pca_res,segment_list,testing_list,breaks_list)
names(final_list)<-c("PCA","Segments","Testing","Breaks")



##### WRITING AS AN OUTPUT THE SEGMENTATION PRODUCTS #####
segment_to_be_written<-data.frame()
for(i in 1:length(segment_list)){
  segment_to_be_written<-rbind.data.frame(segment_to_be_written,segment_list[[i]])  
  
}

write.table(segment_to_be_written,file=paste("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/",name_string,"_segmentation_all_chrs_0_05_breaks_1kb_windows_blaER.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
#pdf(file="length_of_segments_ES_cells.pdf")
#hist(log10(segment_to_be_written[,3]-segment_to_be_written[,2]),main="Length of segments",xlab="log10 (Length)")
#dev.off()

testing_values_to_be_written<-data.frame()
for(i in 1:length(testing_list)){
  temporar<-cbind.data.frame(testing_list[[i]],rep(chr_list[i],length(testing_list[[i]][,1])))
  colnames(temporar)<-c("Mean_Values","Sd_values","Chromosomes")
  testing_values_to_be_written<-rbind.data.frame(testing_values_to_be_written,temporar)  
}

write.table(testing_values_to_be_written,file=paste("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/",name_string,"_mean_and_sd_values_0_05_segmentation_0_05_breaks_blaER.txt",sep=""),quote=F,row.names=F,sep="\t")

breaks_to_be_written<-data.frame()
for(i in 1:length(breaks_list)){
  
  temporar<-cbind.data.frame(breaks_list[[i]],rep(chr_list[i],length(breaks_list[[i]])))
  colnames(temporar)<-c("Break_positions","Chromosomes")
  breaks_to_be_written<-rbind.data.frame(breaks_to_be_written,temporar) 
  
}

write.table(breaks_to_be_written,file=paste("/users/tgraf/aklonizakis/blaer_cells_segmentation/segmentation_results/",name_string,"_breaks_0_05_segmentation_blaER.txt"),quote=F,row.names=F,sep="\t")


