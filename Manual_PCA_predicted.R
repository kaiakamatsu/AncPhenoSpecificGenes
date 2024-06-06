library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(fBasics)
library(MASS)
library(biomaRt)


# --- reread processed tables here 
ge_residuals <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/processed_ge_after_covar.txt")
df <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/transformed_phenos.txt")
covar <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/covars.txt")

# --- plot the residuals on PCA
predicted_ge <- function(x1, x2, y){
  reg = lm( y ~  x1*x2)
  return(reg$fitted.values)
}

ge_fitted <- as.data.frame(apply(ge_residuals, 2, function(x) predicted_ge(bmi, ancestry,x)))
genes <- colnames(ge_fitted)
ge_fitted_transposed <- t(ge_fitted)
pca <- prcomp(ge_fitted_transposed, scale = TRUE, rank = 30)
pc_scores <- as.data.frame(pca$x)
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = "black")+ 
  theme_bw() +
  theme_classic(base_size = 14)

# --- converting common gene names to ENSG
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_list = rownames(sc_gene)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes_list,mart= mart)
genes <- data.frame(genes_list)
colnames(genes) <- c("ensembl_gene_id")
genes_table <- merge(genes, G_list, by = "ensembl_gene_id", all = TRUE)


# --- previously found differential expressed genes for BMI 
de <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/mesa_de.csv")
de <- filter(de, bacon_p < 1.97*10**-6) # threshold set by paper
de_ensembl <- filter(genes_table, hgnc_symbol %in% de$MarkerName)

# --- plot on PCA plots
pc_scores$de <- rep("NO", nrow(pc_scores))
pc_scores$de[which(rownames(pc_scores) %in% de_ensembl$ensembl_gene_id)] <- "YES"
pc_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(data = . %>% filter(de == "NO"), 
             aes(color = de, size = factor(de), alpha = 0.5)) +
  geom_point(data = . %>% filter(de == "YES"), 
             aes(color = de, size = factor(de), alpha = 1)) +
  scale_size_manual(values = c("YES" = 2, "NO" = 1)) +
  scale_alpha_identity() +
  theme_bw() +
  theme_classic(base_size = 14)

ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/PC_de.png", width = 5, height = 5)


# --- TWAS BMI genes 
twas <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/twas.dat")
unique(twas$PANEL) # a tone of tissue panels 
# these expression values are in Monocytes, so work with blood related panels? 
twas_blood <- filter(twas, PANEL == "GTEx.Whole_Blood")
w <- which(twas_blood$TWAS.P < (0.05/nrow(twas_blood))) # bonferroni
twas_blood_significant <- twas_blood[w,]
genes_sig <- sapply(twas_blood_significant$FILE, function(x) strsplit(x, split = "[.]")[[1]][5])
pc_scores$twas <- rep("NO", nrow(pc_scores))
pc_scores$twas[which(rownames(pc_scores) %in% genes_sig)] <- "YES"

pc_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(data = . %>% filter(twas == "NO"), 
             aes(color = twas, size = factor(twas), alpha = 0.5)) +
  geom_point(data = . %>% filter(twas == "YES"), 
             aes(color = twas, size = factor(twas), alpha = 1)) +
  scale_size_manual(values = c("YES" = 2, "NO" = 1)) +
  scale_alpha_identity() +
  theme_bw() +
  theme_classic(base_size = 14)

ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/PC_twas.png", width = 5, height = 5)


