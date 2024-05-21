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
- BMI data: boxcox transformed (lambda = -0.5858586)

## Analysis

### Which genes are differentially expressed in individuals with differnet BMI in different ancestry groups?

> Redundancy analysis

### Comparison to previous association studies

> Cluster based on TWAS results: http://twas-hub.org/traits/

> Cluster based on differential expression analysis results: https://www.nature.com/articles/s41598-019-43881-5
