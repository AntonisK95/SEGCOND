A guide on using SEGCOND
================
2022-09-30

# Introduction

SEGCOND is a custom pipeline that aims to propose genomic regions that
participate into the formation of transcriptional condensates, termed
**P**utative **T**ranscriptional **C**ondensates (PTCs), given a set of
next generation sequencing datasets. The pipeline is flexible and can be
adjusted to your input data.

# Installation and dependencies

Currently SEGCOND is implemented in R as a series of functions. To load
the functions into R, download the file SEGCOND.R in your working
directory and use:

``` r
source("SEGCOND.R")
```

In order to run SEGCOND a series of packages need to be installed first.
A list of them is found below:

1.  strucchange
2.  fitdistrplus
3.  VGAM
4.  gamlss
5.  misha
6.  shaman
7.  Matrix
8.  igraph

Some functions of the algorithm don’t need any of the above packages and
can be run without installing them. An error is going to pop up
specifying the missing package if you try to run a function that needs
it.

For the explanation of every functions’ parameters please refer to the
**Vignette** file.

# Input data files

SEGCOND requires a minimum of two files in order to operate.

1.  A **marks** file containing sorted, continuous genomic bins of fixed
    width and a series of values associated with each bin. A typical
    file would contain average RPKM values of ATAC-seq / ChIP-seq
    experiments per bin. Such files can be easily created when
    ChIP-seq/ATAC-seq data is available using various tools such as
    [deeptools.](https://deeptools.readthedocs.io/en/develop/)

2.  A shaman **track** file. For details on how to obtain such a track
    file from a Hi-C dataset, please refer to the documentation of
    [**shaman**](https://tanaylab.bitbucket.io/shaman/index.html). Note
    that a **track** file is only needed for the last step of SEGCOND.

# Breakdown of the algorithm and toy data example

The algorithm operates on three main steps. In the **data** folder you
can find a toy **marks** file and a **track** file that can be used in
order to run the algorithm.

For more details on the theoretical background and the rationale behind
each step, refer to our publication “**SEGCOND predicts putative
transcriptional condensate- associated genomic regions by integrating
multi-omics data**”.

## Segmentation

The first step of the pipeline separates a given chromosome into
**segments** by employing a time series data analysis algorithm. The
algorithm identifies breakpoints along chromosomes, seperating regions
into segments. Each segment is expected to exhibit different properties
(i.e. different mark values) compared to its neighboring segments.

To run this step we will use the **segment_genome** function on a
dataset provided.

First let’s read the dataset:

``` r
toy_data<-read.delim(sep="\t",header=F,file="chr1_pca_values_2000_bins.txt")

print(toy_data[1:6,])
```

    ##     V1    V2    V3         V4
    ## 1 chr1     1  5000 -2.3304102
    ## 2 chr1  5001 10000 -2.1579603
    ## 3 chr1 10001 15000  1.3037624
    ## 4 chr1 15001 20000  0.9198994
    ## 5 chr1 20001 25000  2.1143079
    ## 6 chr1 25001 30000 12.4207187

The file consists of 2000 continuous 5kb bins spanning from the start of
chr1. Each bin is associated with a PCA value. The PCA value is the
first principal component value (PC1), after applying PCA on a set of
RPKM values associated with each bin. The original file before the PCA
transformation can be found in the **source** files folder and is going
to be used in a later step.

**Important Notes** 

1. The segmentation step process one-dimensional data and that is why a PCA transformation is applied. Since only the first principal component values (PC1) are used, the amount of variance captured by PCA has to be evaluated by the user. In our study, where transcription activation-related marks were used, the PC1 captured more than 50% of the variance. We would recommend using a dataset that exceeds this amount of variance if activation-related marks are used. 

2. Users should also evaluate whether different input datasets contribute equally towards the PC1 values. If they don't, a specfic dataset could be "masking" others. Such practical examples could be two different transcription factor ChIP-seq datasets, with one dominating PC1 contribution while the other not. In such cases it is advised to run the segmentation step for each dataset separately.   

We can partition this part of chr1 into segments by running the
following command:

``` r
segments<-segment_genome(data=toy_data,chr_list="chr1",window=1000,minimum_percentage=0.05,break_threshold=5,run_pca=F)
```

    ## Loading required package: strucchange

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## Loading required package: sandwich

    ## [1] "Now segmenting chr1 !"

``` r
print(segments[[1]])
```

    ##                       chrom   temp1    temp2
    ## chr1 1 625000          chr1       1   625000
    ## chr1 625001 2590000    chr1  625001  2590000
    ## chr1 2590001 3510000   chr1 2590001  3510000
    ## chr1 3510001 3915000   chr1 3510001  3915000
    ## chr1 3915001 5615000   chr1 3915001  5615000
    ## chr1 5615001 6635000   chr1 5615001  6635000
    ## chr1 6635001 6790000   chr1 6635001  6790000
    ## chr1 6790001 8700000   chr1 6790001  8700000
    ## chr1 8700001 10000000  chr1 8700001 10000000

Note that this function may require a lot of processing time depending
on the window and the minimum_percentage options. Smaller windows and a
bigger minimum_percentage will lead to less processing time : However
the decision over this parameters is not trivial and should be performed
with caution.

The window option controls the size of a sliding window which is
employed during the breakpoints generation. This in turn determines the
number of values that are evaluated for each iteration of the algorithm.
Smaller window values will lead to the inspection of less values making
the breakpoint decision more unstable. On the contrary, bigger window
values will lead to unfeasible computing times.

The minimum_percentage defines the minimum length of segments, expressed
as a percentage of the window size, which affects a lot the final
result. Keep in mind the question of interest you are trying to answer
when defining a minimum_percentage cuttoff.

The variable **segments** contains a list with the following slots:

1.  The input data used for the segmentation pipeline. Data is provided
    after the PCA transformation.
2.  Segment coordinates per chromosome.
3.  Average values and standard deviation of the PC1 values within each
    segment.
4.  A set with the regions where a structural break was found and
    reported. This vector can be useful in case you want to plot your
    results.
5.  Segment coordinates merged together as a single data.frame. This
    object is going to be used in downstream steps.

## Annotation

In our case, after segmenting the genome we want to identify segments
that show an abundance of activation-related marks. To do so, we will
need:

1.  A **marks** file. We are going to use the **marks** file that we
    used in the previous step before the PCA transformation.
2.  The segments produced by the previous function.

To run the annotation step we need a file indicating which of the bins
that we used in the segmentation step are actual enhancer elements.This
file can be provided by the user or be prepared from the **marks** file.
Here we are going to prepare it with the **marks** file.

The below function z-transforms each column containing values after log2
transforming them.

``` r
toy_marks<-read.delim(sep="\t",header=T,stringsAsFactors = F,file="marks_0h.txt")

z_score_normalized<-get_z_score(toy_marks)
print(z_score_normalized[1:6,])
```

    ##   X..chr. X.start. X.end. X.ATAC0_1. X.ATAC0_2. X.CEBP0.1. X.CEBP0.2.
    ## 1    chr1        1   5000  -2.967375  -2.992011 -3.4128994  -3.408674
    ## 2    chr1     5001  10000  -2.108887  -1.702335 -2.5799643  -1.836085
    ## 3    chr1    10001  15000  -1.070412  -1.998669 -0.6704501  -1.174664
    ## 4    chr1    15001  20000  -2.225256  -1.476608 -1.2789003  -1.574022
    ## 5    chr1    20001  25000  -2.967375  -1.845077 -2.1968511  -2.825453
    ## 6    chr1    25001  30000  -1.587042  -1.568002 -3.4128994  -3.408674
    ##   X.H000H3K27acX1.pileup_signal. X.H000H3K27acX2.pileup_signal.
    ## 1                      -2.907909                      -2.584390
    ## 2                      -1.885490                      -1.729614
    ## 3                       1.766255                       1.929492
    ## 4                       1.963535                       2.178137
    ## 5                       2.490124                       2.494125
    ## 6                       3.407897                       3.241702
    ##   X.H000H3K4me3X1.pileup_signal. X.H000H3K4me3X2.pileup_signal.
    ## 1                      -2.387807                      -2.593350
    ## 2                      -1.546090                      -1.464726
    ## 3                       2.296942                       2.314732
    ## 4                       1.894107                       1.212814
    ## 5                       2.372626                       1.852075
    ## 6                       4.712924                       4.833217

Using this z-score normalized data.frame, we can come up with genomic
bins that display an abundance of activatory-related marks. Here we deem
any bin that shows a high Z-score for both H3K27ac ChIP-seq signal and
ATAC-seq signal as activation-related.

To get these bins we can run:

``` r
enhancer_bins<-votes_matrix(cooccupancy_mat=z_score_normalized,picked_columns=c(4,5,8,9),threshold=1,picked_lines=FALSE,lines=NULL)
print(enhancer_bins[1:6,])
```

    ##    chr start   end co-occupancy
    ## 1 chr1     1  5000            0
    ## 2 chr1  5001 10000            0
    ## 3 chr1 10001 15000            0
    ## 4 chr1 15001 20000            0
    ## 5 chr1 20001 25000            0
    ## 6 chr1 25001 30000            0

We can now identify segments that show an abundance of the identified
enhancer bins. To do so, we will assume that the distribution of
enhancer bins within a segment of size **N** follows a zero-inflated
negative binomial distribution. Given this assumption, we can model a
background distribution for each segment size and evaluate whether there
is an enrichment of enhancer bins or not for every segment of that size.

``` r
fit<-suppressWarnings(fit_zinegbin_model(segments=segments[[1]],new_scores=enhancer_bins,permutations=1000))
```

    ## Loading required package: fitdistrplus

    ## Loading required package: MASS

    ## Loading required package: survival

    ## Loading required package: VGAM

    ## Loading required package: stats4

    ## Loading required package: splines

    ## Loading required package: gamlss

    ## Loading required package: gamlss.data

    ## 
    ## Attaching package: 'gamlss.data'

    ## The following object is masked from 'package:datasets':
    ## 
    ##     sleep

    ## Loading required package: gamlss.dist

    ## Loading required package: nlme

    ## Loading required package: parallel

    ##  **********   GAMLSS Version 5.4-3  **********

    ## For more on GAMLSS look at https://www.gamlss.com/

    ## Type gamlssNews() to see new features/changes/bug fixes.

    ## [1] "Generating random distributions!"
    ## [1] "Calculating actual hits!"
    ## [1] "Fitting Zero-inflated NB models!"

``` r
print(fit[1:6,])
```

    ##                       chr   start     end     log2FC     p-value
    ## chr1 1 625000        chr1       1  625000 -8.7287981 0.682034558
    ## chr1 625001 2590000  chr1  625001 2590000  2.2106220 0.003433315
    ## chr1 2590001 3510000 chr1 2590001 3510000 -0.5272157 0.416570657
    ## chr1 3510001 3915000 chr1 3510001 3915000  1.5669767 0.032283112
    ## chr1 3915001 5615000 chr1 3915001 5615000 -9.9544869 0.857952684
    ## chr1 5615001 6635000 chr1 5615001 6635000  1.3747649 0.065295268

Note that ocassionally the function may fail to run with an error of
**“mle failed to estimate the parameters with the error code 7”**. In
that case during the model fitting process the parameters of a model
were not calculated. Running the command again or increasing the
permutation parameter usually resolves this issue.

Each segment is now associated with a p-value and an enrichment score.
P-values evaluate how extreme the observed number of enhancer bins found
in each segment are and logFC values show the observed/expected ratio of
enhancer bins for each segment.

## Hi-C Integration

As of now we have partitioned part of chr1 into segments and pinpoint
those that show an abundance of enhancer elements. We will now include
Hi-C data into our analysis, in order to identify segment hubs that
interact frequently in 3D space.

This step needs a shaman **track** in order to be run. The track
consists of point to point normalized interaction scores between
chromosomal coordinates and represents a bin-free normalization
approach. By pooling all interaction scores between a pair of segments,
we can calculate the median/mean value of those interactions.

Note that this procedure is applied for single segments as well. Since
segments are quite large, multiple 3D contacts can occur within the
regions found in the segment.

To run the below function we need the segments that we produced in the
above steps. To minimize computational time, only segments that are
enriched in enhancer bins are going to be passed to the function.

``` r
hic_scores<-get_segment_interactions(misha_path="hg38",hic_normalized_track="blaer_normalized_0h",chr_list="chr1",segments=fit[which(fit[,5]<=0.05),c(1,2,3)],distance_filter=2000000)
```

The resulting list of the above command was saved as an RDS object. We
will now read the RDS object:

``` r
hic_scores<-readRDS("shaman_hic_scores.RDS")
hic_scores
```

    ## $`Mean Scores`
    ## $`Mean Scores`$chr1
    ##                      chr1 625001 2590000 chr1 3510001 3915000
    ## chr1 625001 2590000             8.886967             17.76436
    ## chr1 3510001 3915000           17.764357             14.72115
    ## 
    ## 
    ## $`Median Scores`
    ## $`Median Scores`$chr1
    ##                      chr1 625001 2590000 chr1 3510001 3915000
    ## chr1 625001 2590000                 16.0                 22.6
    ## chr1 3510001 3915000                22.6                 19.5

The hic_scores list now contains two slots. Each slot contains a
symmetric matrix with the segment names as row-names and column-names.
Each entry of the metrix is either the mean (slot 1) or the median (slot
2) score of the normalized interactions between two segments.

## Final identification of PTCs

Using the above information we can isolate segment hubs.

The idea here is to treat each segment as a node in a network. Edges of
the network are only formed if a segment interacts highly in 3D space
with another segment or itself, as determined by the median/mean SHAMAN
score. Moreover each node can be “active” or “inactive” based on the
enhancer content of it. We are interested in:

1.  Segments that are “active”.
2.  Segments that show high SHAMAN scores.

Using the below function, we can isolate segment hubs that exhibit these
properties. Notice again how a hub can consist of a single segment: Such
single segments are very abundant in enhancers and show high normalized
interaction scores over their range, suggesting high enhancer activity.

``` r
final_ptcs<-get_putative_condensates(shaman_hic_scores=hic_scores[[2]],shaman_threshold=17,qualifying_segments=rownames(fit)[which(fit[,5]<=0.05)],distance_filter=2000000)
```

    ## Loading required package: igraph

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
head(final_ptcs[[1]])
```

    ## [[1]]
    ## [1] "chr1 625001 2590000"  "chr1 3510001 3915000"

``` r
head(final_ptcs[[2]])
```

    ## [[1]]
    ##                      chr1 625001 2590000 chr1 3510001 3915000
    ## chr1 625001 2590000                    0                    1
    ## chr1 3510001 3915000                   1                    1

final_ptcs contain a list with two slots. The first slot contains the
identified PTCs. The second slot contains the network node/edge
information for hubs that consist of \>=2 segments.
