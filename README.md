# AncPhenoSpecificGenes

Final Project for BGGN 273, Biostatistics

## Dependencies

R packages

- vegan
- data.table
- dplyr
- ggplot2
- fBasics
- MASS
- biomaRt

## Data Processing

### DATA

Data is extracted from the Multi-Ethnic Study of Atherosclerosis (MESA) study. Individuals are subset to those with both gene expression and phenotype data.

### COVARS

I regress out several confounders/covariates from the gene expression data before continuing with the analysis

- B cell, T cell, NK cell, Neutrophil contamination
- gender
- site of data collection

### DATA TRANSFORMATION

- Gene Expression data: inverse normal transformed
- BMI data: boxcox transformed (lambda = -0.5858586), Z-transformed before RDA analysis

## Analysis

### Which genes are differentially expressed in individuals with differenet BMI? Which genes are differentially expressed in different BMI groups of different ancestries? 

> Redundancy analysis

### What are the functions of these genes? 

### Interesting Findings 

1. ENSG00000135842, NIBAN

> regulates adipocyte apoptosis, important in obesity development (https://pubmed.ncbi.nlm.nih.gov/33185830/)



### Comparison to previous association studies

> Cluster based on TWAS results: http://twas-hub.org/traits/

> Cluster based on differential expression analysis results:
- https://www.nature.com/articles/s41598-019-43881-5 (see table 2)
- https://www.nature.com/articles/s41366-022-01240-x#Sec8 (supplmental table 2)
