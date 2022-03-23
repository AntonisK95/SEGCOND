# Prediction of genomic regions participating in the formation of transcriptional condensates

This repository contains scripts related to the 3 main steps of a computational pipeline that aims to propose genomic regions that participate in the formation of transcriptional condensates. 

More details for this pipeline may be found in our paper:  
**"Identification of transcriptional condensate-associated genomic regions through the integration of multi-omics data"**  
_Antonis Klonizakis, Christoforos Nikolaou and Thomas Graf_

Our methodology aims to identify potential candidate regions participating in the formation and maintenance of biomolecular condensates, using a series of omics experimental outputs, coupled with conformational data such as HiC.
 
Regions in the genome participating in putative, transcriptional condensates should exhibit specific characteristics within regards to their protein occupancy, accessibility levels and 3D chromatin contacts. There are several observations being reported or theorized for regions identified as transcriptional condensate. Condensate regions are reported to be co-occupied by multiple transcription factors and high levels of enhancer related marks, while also harboring multiple 3D interactions between gene promoters and regulatory elements.

The proposed methodology consists of three distinct stages and is outlined below:
Hi-C Integration and candidate identification:  (Figure 1E).


The pipeline can be broken down into 3 stages: 

1. [**Segmentation of the genome**](https://github.com/AntonisK95/Prediction_of_transcr_condensates/tree/main/Segmentation)
   Omics-tracks Integration and Genome Segmentation: We integrate multiple omics datasets and through dimensionality reduction and genome segmentation we create a set of distinct genomic segments in linear chromosomes.

    ![Genome Segmentation](Figures/Figure1.png)

2. [**Annotation of segments**](https://github.com/AntonisK95/Prediction_of_transcr_condensates/tree/main/Annotation)
   Segment Annotation: Each segment is scored and assigned to a different functional class with the focus being on enhancer-associated properties.

   ![Genome Segment Annotation](Figures/Figure2.png)

3. [**Integration of Hi-C data**](https://github.com/AntonisK95/Prediction_of_transcr_condensates/tree/main/Hi-C_Integration)
   3D interaction between and within segments is scored with the integration of Hi-C data. Candidate regions are identified through the application of a set of thresholds associated with chromosomal interaction values. This step makes use of [SHAMAN](https://github.com/tanaylab/shaman).

    ![Hi-C Integration](Figures/Figure3.png)

The scripts that perform the above tasks are written in R. In some cases computation was carried out in a computer cluster that uses the Univa GRID engine as batch system. In those cases, an accompanied script is also provided. 

