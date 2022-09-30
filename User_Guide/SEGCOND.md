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
**viginnete** file.

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
each step, refer to our publication “**SEGCOND proposes putative
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

head(toy_data)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["V1"],"name":[1],"type":["fct"],"align":["left"]},{"label":["V2"],"name":[2],"type":["int"],"align":["right"]},{"label":["V3"],"name":[3],"type":["int"],"align":["right"]},{"label":["V4"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"chr1","2":"1","3":"5000","4":"-2.3304102","_rn_":"1"},{"1":"chr1","2":"5001","3":"10000","4":"-2.1579603","_rn_":"2"},{"1":"chr1","2":"10001","3":"15000","4":"1.3037624","_rn_":"3"},{"1":"chr1","2":"15001","3":"20000","4":"0.9198994","_rn_":"4"},{"1":"chr1","2":"20001","3":"25000","4":"2.1143079","_rn_":"5"},{"1":"chr1","2":"25001","3":"30000","4":"12.4207187","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

The file consists of 2000 continuous 5kb bins spanning from the start of
chr1. Each bin is associated with a PCA value. The PCA value is the
first principal component value (PC1), after applying PCA on a set of
RPKM values associated with each bin. The original file before the PCA
transformation can be found in the **source** files folder and is going
to be used in a later step.

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
head(segments[[1]])
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["chrom"],"name":[1],"type":["fct"],"align":["left"]},{"label":["temp1"],"name":[2],"type":["int"],"align":["right"]},{"label":["temp2"],"name":[3],"type":["int"],"align":["right"]}],"data":[{"1":"chr1","2":"1","3":"625000","_rn_":"chr1 1 625000"},{"1":"chr1","2":"625001","3":"2590000","_rn_":"chr1 625001 2590000"},{"1":"chr1","2":"2590001","3":"3510000","_rn_":"chr1 2590001 3510000"},{"1":"chr1","2":"3510001","3":"3915000","_rn_":"chr1 3510001 3915000"},{"1":"chr1","2":"3915001","3":"5615000","_rn_":"chr1 3915001 5615000"},{"1":"chr1","2":"5615001","3":"6635000","_rn_":"chr1 5615001 6635000"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

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
head(z_score_normalized)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["X..chr."],"name":[1],"type":["chr"],"align":["left"]},{"label":["X.start."],"name":[2],"type":["int"],"align":["right"]},{"label":["X.end."],"name":[3],"type":["int"],"align":["right"]},{"label":["X.ATAC0_1."],"name":[4],"type":["dbl"],"align":["right"]},{"label":["X.ATAC0_2."],"name":[5],"type":["dbl"],"align":["right"]},{"label":["X.CEBP0.1."],"name":[6],"type":["dbl"],"align":["right"]},{"label":["X.CEBP0.2."],"name":[7],"type":["dbl"],"align":["right"]},{"label":["X.H000H3K27acX1.pileup_signal."],"name":[8],"type":["dbl"],"align":["right"]},{"label":["X.H000H3K27acX2.pileup_signal."],"name":[9],"type":["dbl"],"align":["right"]},{"label":["X.H000H3K4me3X1.pileup_signal."],"name":[10],"type":["dbl"],"align":["right"]},{"label":["X.H000H3K4me3X2.pileup_signal."],"name":[11],"type":["dbl"],"align":["right"]}],"data":[{"1":"chr1","2":"1","3":"5000","4":"-2.967375","5":"-2.992011","6":"-3.4128994","7":"-3.408674","8":"-2.907909","9":"-2.584390","10":"-2.387807","11":"-2.593350","_rn_":"1"},{"1":"chr1","2":"5001","3":"10000","4":"-2.108887","5":"-1.702335","6":"-2.5799643","7":"-1.836085","8":"-1.885490","9":"-1.729614","10":"-1.546090","11":"-1.464726","_rn_":"2"},{"1":"chr1","2":"10001","3":"15000","4":"-1.070412","5":"-1.998669","6":"-0.6704501","7":"-1.174664","8":"1.766255","9":"1.929492","10":"2.296942","11":"2.314732","_rn_":"3"},{"1":"chr1","2":"15001","3":"20000","4":"-2.225256","5":"-1.476608","6":"-1.2789003","7":"-1.574022","8":"1.963535","9":"2.178137","10":"1.894107","11":"1.212814","_rn_":"4"},{"1":"chr1","2":"20001","3":"25000","4":"-2.967375","5":"-1.845077","6":"-2.1968511","7":"-2.825453","8":"2.490124","9":"2.494125","10":"2.372626","11":"1.852075","_rn_":"5"},{"1":"chr1","2":"25001","3":"30000","4":"-1.587042","5":"-1.568002","6":"-3.4128994","7":"-3.408674","8":"3.407897","9":"3.241702","10":"4.712924","11":"4.833217","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

Using this z-score normalized data.frame, we can come up with genomic
bins that display an abundance of activatory-related marks. Here we deem
any bin that shows a high Z-score for both H3K27ac ChIP-seq signal and
ATAC-seq signal as activation-related.

To get these bins we can run:

``` r
enhancer_bins<-votes_matrix(cooccupancy_mat=z_score_normalized,picked_columns=c(4,5,8,9),threshold=1,picked_lines=FALSE,lines=NULL)
head(enhancer_bins)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["chr"],"name":[1],"type":["chr"],"align":["left"]},{"label":["start"],"name":[2],"type":["int"],"align":["right"]},{"label":["end"],"name":[3],"type":["int"],"align":["right"]},{"label":["co-occupancy"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"chr1","2":"1","3":"5000","4":"0","_rn_":"1"},{"1":"chr1","2":"5001","3":"10000","4":"0","_rn_":"2"},{"1":"chr1","2":"10001","3":"15000","4":"0","_rn_":"3"},{"1":"chr1","2":"15001","3":"20000","4":"0","_rn_":"4"},{"1":"chr1","2":"20001","3":"25000","4":"0","_rn_":"5"},{"1":"chr1","2":"25001","3":"30000","4":"0","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

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
head(fit)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["chr"],"name":[1],"type":["fct"],"align":["left"]},{"label":["start"],"name":[2],"type":["int"],"align":["right"]},{"label":["end"],"name":[3],"type":["int"],"align":["right"]},{"label":["log2FC"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p-value"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"chr1","2":"1","3":"625000","4":"-8.5978446","5":"0.678929025","_rn_":"chr1 1 625000"},{"1":"chr1","2":"625001","3":"2590000","4":"2.2116693","5":"0.003148344","_rn_":"chr1 625001 2590000"},{"1":"chr1","2":"2590001","3":"3510000","4":"-0.4789136","5":"0.395388327","_rn_":"chr1 2590001 3510000"},{"1":"chr1","2":"3510001","3":"3915000","4":"1.6903028","5":"0.032863423","_rn_":"chr1 3510001 3915000"},{"1":"chr1","2":"3915001","3":"5615000","4":"-9.9283850","5":"0.866989427","_rn_":"chr1 3915001 5615000"},{"1":"chr1","2":"5615001","3":"6635000","4":"1.4144904","5":"0.050854077","_rn_":"chr1 5615001 6635000"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

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
