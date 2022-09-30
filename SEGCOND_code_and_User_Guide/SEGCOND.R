                                      ################################
                                      #### SEGMENTATION ALGORITHM ####
                                      ################################

                         ###########################################################
                         #### NEEDED FUNCTIONS FOR THE RUNNING OF THE ALGORITHM ####
                         ###########################################################

### MERGING SIMILAR BREAKS : IF A BREAK IS X POSITIONS NEXT TO ANOTHER, THE SMALLEST IS TO BE KEPT ###
clearence<-function(input_breaks,break_threshold){
  j<-1
  
  backup<-input_breaks
  while(j<=length(backup)-1){
    
    check<-backup[j+1]-backup[j]
    if(check>=break_threshold){
      j<-j+1  
      
    }
    else{
      
      backup<-backup[-(j+1)]
    }
  }
  return(backup)
}

#### DRAW FINAL COORDINATES OF SEGMENTS BASED ON THE GENERATED BREAKS ####

merge_segments<-function(choosen,chrom,data1,size){
  
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

##### THIS FUNCTION CALCLUATES THE MEAN AND SD OF PCA/PROVIDED VALUES WITHIN EACH SEGMENT ###
variability_domains<-function(domains,data1){
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

                       #################################################################
                       #### MAIN FUNCTION : CALCULATE BREAKPOINTS ACROSS THE GENOME ####
                       #################################################################

segment_genome<-function(data,chr_list,window=1000,minimum_percentage=0.05,break_threshold=5,run_pca=T){
  
### Load strucchange. Return an error message if package is not available. ###

if(require(strucchange)){
  library(strucchange)
}else{
  return(print("strucchange is not installed, please install it first!"))
}

#### LOAD DATA . Format (tab separated): chr start end value-1 ... value-n #####

  
if(run_pca==T){  
  
### Remove "NA" values from the data matrix. ###

data[is.na(data)]<-0


### REDUCE DIMENSIONS OF DATASETS TO A SINGLE VARIABLE ###

print("Running PCA")

test<-prcomp(data[,4:ncol(data)],scale=T)

### IMPORTANT STEP : CHECK THE AMOUNT OF VARIANCE CAPTURED BY THE FIRST PRINCIPLE COMPONENT ###

pca_res<-summary(test)

print(pca_res)

### PREPARE DATA ###

final<-test$x[,1]
final<-cbind.data.frame(data[,1:3],final)
data<-final
colnames(data)<-c("chrom", "start", "end", "pca1")
}

### LISTS TO BE FILLED WITH RESULTS ###
l<-1
segment_list<-list();testing_list<-list();breaks_list<-list()

### INITIATE MAIN SCRIPT ###

for(ch in chr_list){
  grep(paste("^",ch,"$",sep=""), data[,1])->i;data1<-data[i,c(1,2,3,4)]
  
  print(paste("Now segmenting",ch,"!"))
  
  ### CREATE TIME SERIES OBJECT ###
  tsfc<-ts(data1[,4])
  size<-length(data1[,1])
  
  ##### DECIDE THE SIZE OF A SLIDING WINDOW : TESTED 500, 750, 1000 BINS AND 1000 BINS SEEMED TO OPERATE WELL ####
  
  ends<-seq(window,size,round(window/2))
  starts<-seq(1,size-window,round(window/2))
  final_box<-c(size-window,size)
  breaks<-c()
  
  ### IMPORTANT : THE h PARAMETER IS THE MINIMUM SIZE A SEGMENT SHOULD HAVE. THIS IS A NON TRIVIAL PARAMETER ###
  ### THAT HAS TO BE EVALUATED SOMEHOW. AFTER TESTING, A 5% VALUE SEEMED APPROPRIATE ####
  
  for(i in 1:length(starts)){
    testing<-tsfc[starts[i]:ends[i]]
    breakpoints(testing ~ 1, h=minimum_percentage)->bp;
    results<-bp$breakpoints
    results<-results+(starts[i]-1)
    breaks<-c(breaks,results)
    
  }
  testing<-tsfc[final_box[1]:final_box[2]]
  breakpoints(testing ~ 1, h=minimum_percentage)->bp;
  results<-bp$breakpoints
  results<-results+(final_box[1]-1)
  breaks<-c(breaks,results)
  
  #### THE BREAKS VECTOR CONTAINS NON FILTERED BORDERS CALCULATED BY THE SLIDING WINDOW APPROACH ####
  
  breaks<-sort(unique(breaks))
  
  #### MERGE BREAKS THAT ARE VERY CLOSE TO EACH OTHER / IDENTICAL 
  breaks<-clearence(breaks,break_threshold = break_threshold)
  
  #### FILLING LISTS WITH RESULTS ! ####
  
  final<-merge_segments(breaks,chrom=ch,data1=data1,size=size)
  testing<-variability_domains(final,data1=data1)
  segment_list[[l]]<-final
  testing_list[[l]]<-testing
  breaks_list[[l]]<-breaks
  l<-l+1
  
}


names(segment_list)<-chr_list
names(testing_list)<-chr_list
names(breaks_list)<-chr_list

for(k in 1:length(segment_list)){
  
  if(k==1){
    
    final_frame<-segment_list[[k]]
    
  }
  
  else{
    
    final_frame<-rbind.data.frame(final_frame,segment_list[[k]])
  }
  
}
rownames(final_frame)<-paste(final_frame[,1],final_frame[,2],final_frame[,3])
final_list<-list(final_frame,segment_list,testing_list,breaks_list,data)
names(final_list)<-c("All Segments","Segments per Chromosome","Mean_SD_values_of_segments","Breaks","Data used")

return(final_list)

}

                      #################################################################
                      ############### SEGMENT ANNOTATION FUNCTIONS ####################
                      #################################################################


                  ######################################################################
                  #### NEEDED FUNCTIONS FOR THE RUNNING OF THE ANNOTATION ALGORITHM ####
                  ######################################################################



###### DATA PREPARATION ######
### In case you have a file ready you can skip these functions  ####
### Necessary file : chromosome - start - end - 0/1 column ####
### The chromosomal coordinates have to be the same with the ones used to carry out the segmentation. ###


#### File Format : chr\tstart\tend\tmark1\tmark2....markN ####
#### THIS FUNCTION LOG2 TRANSFORMS COLUMNS AFTER ADDING A 0.01 PSEDUCOUNT ####
#### AFTERWARDS IT Z-TRANSFORMS THE COLUMNS ####

get_z_score<-function(marks){
  backup<-marks
  backup[is.na(backup)]<-0

  for(j in 4:ncol(backup)){
    
    means_of_log<-mean(log2(backup[,j]+0.01))
    sds_of_log<-sd(log2(backup[,j]+0.01))
    backup[,j]<-(log2(backup[,j]+0.01)-means_of_log)/sds_of_log
  }
  
  return(backup)
  
}

### Transforms z-scaled columns into "interesting" bins : Every bin that has a Z score of >=threshold for a set amount of marks ###

votes_matrix<-function(cooccupancy_mat,picked_columns=NULL,threshold=NULL,picked_lines=FALSE,lines=NULL){
  
  if(picked_lines==FALSE){
    lines_of_interest<-as.numeric(which(unlist(apply(cooccupancy_mat[,c(picked_columns)]>=threshold,1,function(x){if(length(grep(pattern=FALSE,x))==0){return(TRUE)}else{return(FALSE)}}))))
  }
  else{
    
    lines_of_interest<-lines
  }
  new_test<-cbind.data.frame(cooccupancy_mat[,c(1,2,3)],rep(0,length(cooccupancy_mat[,1])))
  colnames(new_test)<-c("chr","start","end","co-occupancy")
  new_test[lines_of_interest,4]<-1
  return(new_test)
}


           ############################################################################################################
           #### MAIN FUNCTIONS : 1. RANDOM DISTRIBUTION, 2. ACTUAL DISTRIBUTION, 3. FITTING ZERO-INFLATED NB MODEL ####
           ############################################################################################################

randomize_segment_votes<-function(segments,new_scores,permutations){
  
  mean_random_votes_per_segment_cooccupancy<-list()
  permutation_votes<-list()
  ### Calculate size distribution of segments. ###
  
  segment_length<-unique(sort(segments[,3]-segments[,2]))
  binsize<-new_scores[1,3]-new_scores[1,2]+1
  for(j in 1:length(segment_length)){
    test<-segment_length[j]
    test2<-round((test)/binsize)
    
    mean_vector<-vector()
    
    random_starts<-sample(seq(1,length(new_scores[,1])-test2,1),permutations)
    random_ends<-random_starts+test2
    
    for(i in 1:length(random_starts)){
      mean_vector[i]<-sum(new_scores[random_starts[i]:random_ends[i],4])
    }
    permutation_votes[[j]]<-mean_vector
  }
  
  names(permutation_votes)<-segment_length
  
  for(i in 1:nrow(segments)){
    
    segment_size<-segments[i,3]-segments[i,2]
    mean_random_votes_per_segment_cooccupancy[[i]]<-permutation_votes[[as.character(segment_size)]]
  }
  
  
  return(mean_random_votes_per_segment_cooccupancy)
}


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

#### Fit a zero-inflated negative binomial distribution to background distribution ####
#### Get a p-value and an enrichment score for each segment ####

get_zinegbin_pvalues<-function(segments,mean_random_votes_per_segment_cooccupancy,actual_cooccupancy_votes){
  
  
  if(require(fitdistrplus)){
    library(fitdistrplus)
  }else{
    return(print("fitdistrplus is not installed, please install it first!"))
  }
  
  if(require(VGAM)){
    library(VGAM)
  }else{
    return(print("VGAM is not installed, please install it first!"))
  }
  
  if(require(gamlss)){
    library(gamlss)
  }else{
    return(print("gamlss is not installed, please install it first!"))
  }
  
  
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
  
  to_return<-cbind.data.frame(segments,enrichment_cooccupancy,pvalues_cooccupancy);colnames(to_return)<-c("chr","start","end","log2FC","p-value")
  return(to_return)
}

### PENDING : IN CASE A MODEL FIT, SKIP THAT MODEL FIT AND PROCCEED WITH THE NEXT SEGMENT 
### PENDING : IMAGES TO ASSESS QUALITY CONTROL OF FIT FOR GENERATED RANDOM MODELS. 

                                       ##### WRAP-UP FUNCTION ######

fit_zinegbin_model<-function(segments,new_scores,permutations){
  if(require(fitdistrplus)){
    library(fitdistrplus)
  }else{
    return(print("fitdistrplus is not installed, please install it first!"))
  }
  
  if(require(VGAM)){
    library(VGAM)
  }else{
    return(print("VGAM is not installed, please install it first!"))
  }
  
  if(require(gamlss)){
    library(gamlss)
  }else{
    return(print("gamlss is not installed, please install it first!"))
  }
  print("Generating random distributions!")
  random_votes_per_segment_cooccupancy<-randomize_segment_votes(segments = segments,new_scores = new_scores,permutations = permutations)
  
  print("Calculating actual hits!")
  actual_cooccupancy_votes<-get_actual_votes(segments=segments,new_scores=new_scores)
  
  print("Fitting Zero-inflated NB models!")
  final<-get_zinegbin_pvalues(segments = segments,mean_random_votes_per_segment_cooccupancy = random_votes_per_segment_cooccupancy,actual_cooccupancy_votes = actual_cooccupancy_votes)
  rownames(final)<-paste(final[,1],final[,2],final[,3])
  return(final)
}



                                        ################################
                                        ####  HI-C INTEGRATION STEP ####
                                        ################################

#### THE MAIN COMMAND ACCEPTS A SET OF SEGMENTS COORDINATES AND A LIST OF CHROMOSOMES ###
#### REQUIREMENTS : A SHAMAN NORMALIZED HIC TRACK. ####
#### THE FUNCTION WILL CALCULATE ALL THE SHAMAN NORMALIZED SCORES ASSOCIATED WITH ALL SEGMENTS ####
#### ONE CAN FILTER FOR THE FOLLOWING :
#### YOU CAN ONLY USE AS INPUTS INTERESTING SEGMENTS. AFTER APPLYING A DISTANCE THRESHOLD, SEGMENTS FURTHER AWAY ARE NOT CONSIDERED ANYMORE. ####
#### A SCORE OF 0 IS ASSIGNED TO THOSE. ####

get_segment_interactions<-function(misha_path,hic_normalized_track,chr_list,segments,distance_filter=Inf){
  
  if(require(misha)){
    library(misha)
  }else{
    return(print("misha is not installed, please install it first!"))
  }
  
  if(require(shaman)){
    library(shaman)
  }else{
    return(print("shaman is not installed, please install it first!"))
  }
  
  if(require(Matrix)){
    library(Matrix)
  }else{
    return(print("Matrix is not installed, please install it first!"))
  }
  
  gdb.init(misha_path)
  gsetroot(misha_path)
  gdb.reload()
  
  
  list_with_means<-list();list_with_medians<-list()
  for(i in 1:length(chr_list)){
    temporar<-segments[which(segments[,1]==chr_list[i]),]
    final_matrix<-matrix(nrow=length(temporar[,1]),ncol=length(temporar[,1]))
    final_matrix_median<-matrix(nrow=length(temporar[,1]),ncol=length(temporar[,1]))
    colnames(final_matrix)<-paste(paste(temporar[,1],temporar[,2],temporar[,3]))
    rownames(final_matrix)<-paste(paste(temporar[,1],temporar[,2],temporar[,3]))
    colnames(final_matrix_median)<-paste(paste(temporar[,1],temporar[,2],temporar[,3]))
    rownames(final_matrix_median)<-paste(paste(temporar[,1],temporar[,2],temporar[,3]))
    
    for(k in 1:length(final_matrix[,1])){
      for(j in k:length(final_matrix[,1])){
        distance1<-abs(temporar[k,3]-temporar[j,2])
        distance2<-abs(temporar[k,2]-temporar[j,3])
        distance<-min(distance1,distance2)
        if(distance>distance_filter){
          
          final_matrix[k,j]<-0
          final_matrix_median[k,j]<-0
        }
        else{
        point_score<-gextract(hic_normalized_track, gintervals.2d(temporar[k,1], temporar[k,2], temporar[k,3], temporar[j,1], temporar[j,2], temporar[j,3]), colnames="score")[,7]
        
        if(length(point_score>0)){
          final_matrix[k,j]<-mean(point_score)
          final_matrix_median[k,j]<-median(point_score)
        }
        else{
          final_matrix[k,j]<-0
          final_matrix_median[k,j]<-0
          
        }
        
        
      }
    }}
    list_with_means[[i]]<-final_matrix;list_with_medians[[i]]<-final_matrix_median
  }
  names(list_with_means)<-chr_list;names(list_with_medians)<-chr_list
  list_with_means<-lapply(list_with_means,forceSymmetric);list_with_medians<-lapply(list_with_medians,forceSymmetric);
  list_with_means<-lapply(list_with_means,as.matrix);list_with_medians<-lapply(list_with_medians,as.matrix)
  list_to_return<-list(list_with_means,list_with_medians)
  names(list_to_return)<-c("Mean Scores","Median Scores")
  return(list_to_return)
}


                                   ###### FINAL CONDENSATES CALLING STEP ########

### ADD STEP UP TO THE POINT OF THE GENERATION OF THE CONDENSATES ###
### SHOULD I ADD TAGGING OF EACH SEGMENT AS WELL? ###
### RELATION MATRIX FOR EACH CONDENSATE FROM THE GRAPH? ###

### Given a set of interesting segments and a SHAMAN threshold, this changes all entries into either 0 or 1's ###
conversion_to_binary_matrix<-function(testing_list,qualifiers,threshold){
  for(i in 1:ncol(testing_list)){
    
    ### if a column is not part of the qualifying segments, everything will be turned to 0 ###
    column_name<-colnames(testing_list)[i]
    if(column_name%in%qualifiers){
      
      ones<-which(testing_list[,i]>=threshold)  
      
      ### check if rownames are actually S.E. related ###
      new_ones<-vector()
      
      if(length(ones)==0){
        new_ones<-vector()
      }
      else{
        for(j in 1:length(ones)){
          if(rownames(testing_list)[ones[j]]%in%qualifiers){
            new_ones<-append(new_ones,ones[j])
          }
        }}
      
      testing_list[,i]<-0
      testing_list[new_ones,i]<-1
    }
    else{
      testing_list[,i]<-0
    }}
  return(as.data.frame(as.matrix(testing_list)))
}

distance_filter_and_deletion<-function(testing_for_dist,distance){
  
  for(i in 1:ncol(testing_for_dist)){
    
    rows<-which(testing_for_dist[,i]==1)
    
    
    if(length(rows)>0){
      
      column_name<-colnames(testing_for_dist)[i]
      start_column<-as.numeric(unlist(strsplit(column_name," "))[2])
      end_column<-as.numeric(unlist(strsplit(column_name," "))[3])
      
      all_row_names<-unlist(strsplit(rownames(testing_for_dist)[rows]," "))
      all_row_names_start<-as.numeric(all_row_names[seq(2,length(all_row_names)-1,3)])
      all_row_names_end<-as.numeric(all_row_names[seq(3,length(all_row_names),3)])
      
      diff_1<-abs(start_column-all_row_names_end)
      diff_2<-abs(end_column-all_row_names_start)
      distances_vector<-vector()
      for(j in 1:length(diff_1)){
        
        distances_vector[j]<-min(diff_1[j],diff_2[j])
        
      }
      
      testing_for_dist[rows[which(distances_vector>distance)],i]<-0
      
    }
    
    
  }
  
  
  initial_col<-colnames(testing_for_dist);backup_names<-vector()
  ### We are going to take into advantage the following : In case that a column is zero-inflated that means that the same row is also zero-inflated.
  ### We only need one loop to drop non-connected entitites! 
  for(i in initial_col){
    if(length(testing_for_dist)==1){
      checking<-sum(testing_for_dist)
    }
    else{
      checking<-sum(testing_for_dist[,i])
    }
    if(checking==1){
      backup_names<-append(backup_names,i)
    }
    
    if(checking==0){
      ### NOTE that the matrix is symmetric : During the binarization step, I included a duplication event that copies the upper diagonal
      ### of the matrix to the lower diagonal. Thus if a colname is zero-inflated, the row with the same name will follow that trend
      ### in some cases the data frame shrinks completly and stops being a data frame! Have to take into account this 
      ### Notice that if there are 0 lines the whole thing will have to disappear . If one line is surviving, then a single "1" is going to be returned as output! ####
      
      if(length(testing_for_dist)==1){
        testing_for_dist<-data.frame()
        return(testing_for_dist)
      }
      
      
      testing_for_dist<-testing_for_dist[!(names(testing_for_dist) %in% i),!(names(testing_for_dist) %in% i)]
    }
    
  }
  if(length(testing_for_dist)>1){
    return(testing_for_dist)
  }
  
  else{
    names(testing_for_dist)<-backup_names
    return(testing_for_dist)
  }
}

clean_components<-function(toy_example){
  if(require(igraph)){
    library(igraph)
  }else{
    return(print("igraph is not installed, please install it first!"))
  }
  if(length(toy_example)!=0){
    
    testing_graph<-graph_from_adjacency_matrix(as.matrix(toy_example))
    connected_components<-components(testing_graph)
    list_with_comp<-list()
    for(i in 1:connected_components$no){
      
      list_with_comp[[i]]<-names(connected_components$membership)[which(connected_components$membership==i)]
      
    }
    
    
    return(list_with_comp)
  }}


#### The above have to be run for every chromosome. Afterwards results can be pooled. ###
#### The final list returned should be appended with all entries. ####


get_putative_condensates<-function(shaman_hic_scores,shaman_threshold,qualifying_segments,distance_filter){
  list_with_putative_condensates<-list()
  list_with_all_info<-list()
  for(i in 1:length(shaman_hic_scores)){
    
    binary_matrix<-conversion_to_binary_matrix(testing_list=shaman_hic_scores[[i]],qualifiers=qualifying_segments,threshold=shaman_threshold)
    final_matrix<-distance_filter_and_deletion(testing_for_dist=binary_matrix,distance=distance_filter)
    ptcs<-clean_components(final_matrix)
    list_with_putative_condensates<-append(list_with_putative_condensates,ptcs)
    temporar_list<-list()
    for(j in 1:length(ptcs)){
    temporar_list[[j]]<-final_matrix[ptcs[[j]],ptcs[[j]]]
    }
    list_with_all_info<-append(list_with_all_info,temporar_list)
  }
  
  
  list_with_ptcs<-list(list_with_putative_condensates,list_with_all_info)
  names(list_with_ptcs)<-c("PTCs","Edge Info")
  return(list_with_ptcs)
}

