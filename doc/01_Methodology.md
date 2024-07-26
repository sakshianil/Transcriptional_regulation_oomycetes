# Methodology

![Methodology for motif discovery in Pl. halstedii](../figures/pipeline_1.png)

## Data Preprocessing
The raw data is cleaned and normalized to remove any inconsistencies.

## Differential Expression Analysis
Differential expression analysis is performed using DESeq2, identifying genes that are significantly up- or down-regulated.

## Clustering
Gene expression profiles are clustered to identify groups of co-expressed genes.

## Motif Discovery
Motif discovery is performed using MEME to find common regulatory elements in the promoter regions of co-expressed genes.



![Methodology for motif discovery in closely related oomycetes](../figures/pipeline_2.png)


![Methodology for identification of transcription factors in oomycetes](../images/pipeline_3.png)
