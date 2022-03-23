# Prediction of genomic regions participating in the formation of transcriptional condensates

This repository contains scripts related to the 3 main steps of a computational pipeline that aims to propose genomic regions that participate in the formation of transcriptional condensates. 

More details for this pipeline may be found in our paper:  
**"Identification of transcriptional condensate-associated genomic regions through the integration of multi-omics data"**  
_Antonis Klonizakis, Christoforos Nikolaou and Thomas Graf_

Our methodology aims to identify potential candidate regions participating in the formation and maintenance of biomolecular condensates, using a series of omics experimental outputs, coupled with conformational data such as HiC.

The 

The pipeline can be broken down into 3 stages: 

1. [**Segmentation of the genome**](https://github.com/AntonisK95/Prediction_of_transcr_condensates/tree/main/Segmentation)
2. [**Annotation of segments**](https://github.com/AntonisK95/Prediction_of_transcr_condensates/tree/main/Annotation)
3. [**Integration of Hi-C data**](https://github.com/AntonisK95/Prediction_of_transcr_condensates/tree/main/Hi-C_Integration)

The scripts that perform the above tasks are written in R. In some cases computation was carried out in a computer cluster that uses the Univa GRID engine as batch system. In those cases, an accompanied script is also provided. 

