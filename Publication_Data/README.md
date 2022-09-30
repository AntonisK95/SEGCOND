Code to reproduce our published data.
================
2022-09-30

This is a step-by-step guide reproducing the final files used in the
analysis of our paper. For details regarding the methodology used here
please refer to our publication, the **User_Guide** folder and the
**Viginnetttes** file.

## Segmentation

This step operates on a set of **marks** files. These files can be found
in the folder **Marks_Files**. They include genomic bins of fixed width
(5000 bp) accompanied by RPKM values from various ChIP-seq experiments
and an ATAC-seq experiment. The experiments were performed in a system
of transdifferentiation previously described in the lab, at timepoints
0d, 1d and 7d.

In order to segment the chromosomes for every timepoint we can run the
following commands:

``` r
source("SEGCOND.R")

marks_0d<-read.delim(sep="\t",stringsAsFactors=F,header=T,file="All_marks_timepoint_0h.txt")
marks_1d<-read.delim(sep="\t",stringsAsFactors=F,header=T,file="All_marks_timepoint_24h.txt")
marks_7d<-read.delim(sep="\t",stringsAsFactors=F,header=T,file="All_marks_timepoint_7d.txt")


segments_0d<-segment_genome(data=marks_0d,chr_list=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),window=1000,minimum_percentage=0.05,break_threshold=5,run_pca=TRUE)

segments_1d<-segment_genome(data=marks_1d,chr_list=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),window=1000,minimum_percentage=0.05,break_threshold=5,run_pca=TRUE)

segments_7d<-segment_genome(data=marks_7d,chr_list=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),window=1000,minimum_percentage=0.05,break_threshold=5,run_pca=TRUE)
```

## Annotation of segments

We can now isolate segments that exhibit an abundance of enhancer bins.
We will consider any bin that show simultaneously high levels of H3K27ac
ChIP-seq signal and ATAC-seq signal as an “enhancer” bin.

We can calculate an enrichment score and a p-value for every segment
with the following commands:

``` r
z_score_normalized_0d<-get_z_score(marks_0d)
z_score_normalized_1d<-get_z_score(marks_1d)
z_score_normalized_7d<-get_z_score(marks_7d)

enhancer_bins_0d<-votes_matrix(cooccupancy_mat=z_score_normalized_0d,picked_columns=c(4,5,8,9),threshold=1,picked_lines=FALSE,lines=NULL)

enhancer_bins_1d<-votes_matrix(cooccupancy_mat=z_score_normalized_1d,picked_columns=c(4,5,8,9),threshold=1,picked_lines=FALSE,lines=NULL)

enhancer_bins_7d<-votes_matrix(cooccupancy_mat=z_score_normalized_7d,picked_columns=c(4,5,8,9),threshold=1,picked_lines=FALSE,lines=NULL)

fit_0d<-suppressWarnings(fit_zinegbin_model(segments=segments_0d[[1]],new_scores=enhancer_bins_0d,permutations=1000))
fit_1d<-suppressWarnings(fit_zinegbin_model(segments=segments_1d[[1]],new_scores=enhancer_bins_1d,permutations=1000))
fit_7d<-suppressWarnings(fit_zinegbin_model(segments=segments_7d[[1]],new_scores=enhancer_bins_7d,permutations=1000))
```

The final files of the generated segments alongside logFC and p-value
scores for every timepoint are saved in the
**Segmentation_Results_with_zinegbin_models** folder.

## Hi-C integration

We can now calculate SHAMAN interaction scores for every segment pair.
To do so a SHAMAN Hi-C **track** was generated for every timepoint. For
instructions on how to generate such a file, please refer to the guide
of **shaman**\[<https://tanaylab.bitbucket.io/shaman/index.html>\]

To calculate the scores we can run:

``` r
hic_scores_0d<-get_segment_interactions(misha_path="hg38",hic_normalized_track="blaer_normalized_0h",chr_list=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),segments=fit_0d,distance_filter=Inf)

hic_scores_1d<-get_segment_interactions(misha_path="hg38",hic_normalized_track="blaer_normalized_1d",chr_list=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),segments=fit_1d,distance_filter=Inf)

hic_scores_7d<-get_segment_interactions(misha_path="hg38",hic_normalized_track="blaer_normalized_7d",chr_list=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"),segments=fit_7d,distance_filter=Inf)
```

The symmetric matrix containing the median scores of the interactions
between all segment pairs, for every timepoint, is stored in the
**SHAMAN_Median_scores_segment_pairs** folder.

## Isolation of PTCs

To finally isolate PTCs for every timepoint we can run the following
commands:

``` r
final_ptcs_0d<-get_putative_condensates(shaman_hic_scores=hic_scores_0d[[2]],shaman_threshold=17,qualifying_segments=rownames(fit_0d)[which(fit_0d[,5]<=0.05 & fit_0d[,4]>=0.5849625)],distance_filter=2000000)

final_ptcs_1d<-get_putative_condensates(shaman_hic_scores=hic_scores_1d[[2]],shaman_threshold=17,qualifying_segments=rownames(fit_1d)[which(fit_1d[,5]<=0.05 & fit_1d[,4]>=0.5849625)],distance_filter=2000000)

final_ptcs_7d<-get_putative_condensates(shaman_hic_scores=hic_scores_7d[[2]],shaman_threshold=17,qualifying_segments=rownames(fit_7d)[which(fit_7d[,5]<=0.05 & fit_7d[,4]>=0.5849625)],distance_filter=2000000)
```

The final segments that were found to be part of **PTCs** can be found
in the **Final_Condensates** folder. Each segment is accompanied by a
counter, indicating its **PTC** membership : If two segments have the
same counter, then both of them are members of a **PTC**.
