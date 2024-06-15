# AncPhenoSpecificGenes

Final Project for BGGN 273, Biostatistics

see BGGN273.pdf for 1 slide summary

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

### BMI vs ancestry

> ANOVA

### Which genes are differentially expressed in individuals with differenet BMI?

> Redundancy analysis without BMI/ancestry interaction

### Which genes are differentially expressed in different BMI groups of different ancestries?

> Redundancy analysis with BMI/ancestry interaction

### What are the functions of these genes?

### Interesting Findings

1. RRAGD, ENSG00000025039 (p = 1.234555e-18, BMI associated)

> - plays a crucial role in the cellular response to amino acid availability through regulation of the mTORC1 signaling cascade (genecards)
> - previous gene-level associations (through genetic data) found with weight, bmi, type 2 diabetes (Dornbos et al. 2022 Cell Metabolism)
> - causative experimentally for kidney disease (Schlingmann et al. 2021 J Am Soc Nephrol)
> - kidney diseases linked to bmi (Herrington et al. 2017 PLoS One)

### Comparison to previous association studies

> Cluster based on TWAS results: http://twas-hub.org/traits/

> Cluster based on differential expression analysis results:

- https://www.nature.com/articles/s41598-019-43881-5 (see table 2)
- https://www.nature.com/articles/s41366-022-01240-x#Sec8 (supplmental table 2)
