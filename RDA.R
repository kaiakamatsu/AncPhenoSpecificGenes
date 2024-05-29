install.packages('vegan')
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(fBasics)
library(MASS)
library(biomaRt)

# --- Read in data
df <- as.data.frame(fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/pheno_expression_data.txt.gz"))
colnames(df)
df <- na.omit(df)

# --- Basic information about data 
nrow(df) # 1111 individuals 
table(df$race)
  # AA CAU HIS 
  #209 556 346
table(df$gender)
  # FEMALE   MALE 
  #  552    559
summary(df$age1c)
  # mean = 60

# --- Extract covariate information / confounders 
# we want to remove the variance explained by cell type contamination (B cell, T cell NK cell, Neutrophil), gender, site collection. 
covar <- df[,c(3:9)]
df <- df[,-c(3:9)]
# set gender and site to integers
unique(covar$gender)
gender_int <- setNames(seq_along(unique(covar$gender)) - 1, unique(covar$gender))
covar <- covar %>%
  mutate(gender = gender_int[gender])
unique(covar$site)
site_int <- setNames(seq_along(unique(covar$site)) - 1, unique(covar$site))
covar <- covar %>%
  mutate(site = site_int[site])

# --- Extract gene expression data 
ge_col <- grep("ENSG", colnames(df))
ge <- df[,ge_col]
df <- df[,-ge_col]

# --- Visualize the distributions of data / Normality Test
histogram_normalTest <- function(x){
  hist(x)
  normalTest(x, "da")
}

histogram_normalTest(df$height_cm)
histogram_normalTest(df$weight_lbs)
histogram_normalTest(df$`bmi_(kg/m^2)`)
histogram_normalTest(df$`ldl_(mg/dl)`)
histogram_normalTest(df$`hdl_(mg/dl)`)
histogram_normalTest(df$`total_chol_(mg/dl)`)
histogram_normalTest(ge$ENSG00000000003) # one gene expression column as example
# none of these columns are normal

# --- normalize gene expression data 
# for gene expression data, inverse normal transform (instead of trying to transform every gene differently)
inverse_normal_transform <- function(df){
  result <- apply(df, 2, function(r) qnorm((rank(r,na.last="keep")-0.5)/sum(!is.na(r))) )
  return(data.frame(result))
}
# used by GTEx project for gene expression data
# https://yingji15.github.io/2019-08-17-inverse-normal-transformation/
# https://www.nature.com/articles/nature11401#Sec4 supplement

# feed into ivt function, replace in dataframe
ge_transformed <- inverse_normal_transform(ge)

# --- normalize other columns
boxcox_transform <- function(BCX, data){
  boxcoxlambda <- BCX$x[which.max(BCX$y)]
  print(boxcoxlambda)
  boxcoxdata <- (data^(boxcoxlambda) - 1)/boxcoxlambda
  return(boxcoxdata)
}

height = df$height_cm
weight = df$weight_lbs
bmi = df$`bmi_(kg/m^2)`
LDL = df$`ldl_(mg/dl)`
HDL = df$`hdl_(mg/dl)`
Total_Chol = df$`total_chol_(mg/dl)`

BCX <- boxcox(lm(height ~ 1))
histogram_normalTest(log(df$height_cm)) # log transform does not work
df$height_cm <- qnorm((rank(df$height_cm,na.last="keep")-0.5)/sum(!is.na(df$height_cm))) # quantile normalize

BCX <- boxcox(lm(weight ~ 1))
histogram_normalTest(log(df$weight_lbs)) # log transform works
df$weight_lbs <- scale(log(df$weight_lbs))

BCX <- boxcox(lm(bmi ~ 1))
histogram_normalTest(boxcox_transform(BCX,df$`bmi_(kg/m^2)`)) # -0.5858586
df$`bmi_(kg/m^2)` <- scale(boxcox_transform(BCX,df$`bmi_(kg/m^2)`))

BCX <- boxcox(lm(LDL ~ 1))
histogram_normalTest(boxcox_transform(BCX,df$`ldl_(mg/dl)`))
df$`ldl_(mg/dl)` <- scale(boxcox_transform(BCX,df$`ldl_(mg/dl)`)) # 0.5454545

BCX <- boxcox(lm(HDL ~ 1))
histogram_normalTest(boxcox_transform(BCX,df$`hdl_(mg/dl)`))
df$`hdl_(mg/dl)` <- scale(boxcox_transform(BCX,df$`hdl_(mg/dl)`)) # -0.2222222

BCX <- boxcox(lm(Total_Chol ~ 1))
histogram_normalTest(boxcox_transform(BCX,df$`total_chol_(mg/dl)`)) # does not work
df$`total_chol_(mg/dl)` <- qnorm((rank(df$`total_chol_(mg/dl)`,na.last="keep")-0.5)/sum(!is.na(df$`total_chol_(mg/dl)`))) # quantile normalize

# --- regress out covariates from gene expression data 
regress_out_covars <- function(y, covars){
  reg = summary(lm( y ~ as.matrix(covars) ))
  return(scale(reg$resid))
}
ge_residuals <- apply(ge_transformed, 2, function(x) regress_out_covars(x, covar))
ge_residuals <- as.data.frame(ge_residuals)

# --- check point, save tables of gene expression values after regressed covariates, table of covariates used, and explanatory variables used
write.table(ge_residuals, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/processed_ge_after_covar.txt", quote = F, row.names = F, col.names = T, sep = '\t')
write.table(df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/transformed_phenos.txt", quote = F, row.names = F, col.names = T, sep = '\t')
write.table(covar, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/covars.txt", quote = F, row.names = F, col.names = T, sep = '\t')

# --- restart here
ge_residuals <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/processed_ge_after_covar.txt")
df <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/transformed_phenos.txt")
covar <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/covars.txt")

# --- run RDA analysis
colnames(df)
# lets first try using race and bmi as the explanatory variables for gene expression
  # this will hopefully help us find genes which are differentially expressed in different bmi/ancestry groups
ancestry <- df$race
bmi <- df$`bmi_(kg/m^2)`
bmi_anc_rda <- rda(formula = ge_residuals ~ ancestry*bmi)
# constrained proportion: 0.01368
# uncinstrained proportion: 0.9863

# --- anova test against null model 
null_rda <- rda(formula = ge_residuals ~ 1)
test <- anova.cca(null_rda, bmi_anc_rda)
# 999 permutations
#Model 1: ge_residuals ~ 1
#Model 2: ge_residuals ~ ancestry * bmi
#   ResDf ResChiSquare Df ChiSquare     F Pr(>F)    
#1  1110        22475                               
#2  1105        22168  5    307.37 3.0644  0.001

# --- plot RDA results 
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_vegan_scale1.png", width = 800, height = 800)
ordiplot(bmi_anc_rda, scaling = 1, type = "points")
dev.off()
# interpretation is difficult, constrained variance is very small so we are only looking at 1.3% of the total variance in this plot

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_vegan_scale2.png", width = 800, height = 800)
ordiplot(bmi_anc_rda, scaling = 2, type = "points")
dev.off()
  # black dots = people (rows)
  # red points = genes

# --- custom RDA plot 
perc <- round(100*(summary(bmi_anc_rda)$cont$importance[2, 1:2]), 2)
sc_ind <- scores(bmi_anc_rda, display="sites", choices=c(1,2), scaling=1)
sc_gene <- scores(bmi_anc_rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(bmi_anc_rda, display="bp", choices=c(1, 2), scaling=1)
  # AA gone??? AA used as reference 

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1.png", width = 800, height = 800)
plot(bmi_anc_rda,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlim = c(-3,3), 
     ylim = c(-3,3),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)"), 
     cex.lab = 1.5,
     cex.axis = 1.2
)
color_ind <- rgb(70/255, 130/255, 180/255, alpha = 0.5)  # steelblue with alpha
color_gene <- rgb(242/255, 189/255, 51/255, alpha = 0.5)
points(sc_ind, 
       pch = 21,
       col = "white",
       bg = color_ind,
       cex = 1)
points(sc_gene, 
       pch = 22,
       col = "white",
       bg = color_gene, 
       cex = 1)
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*40, sc_bp[,2]*40, # scale up the score values
       col = "red", 
       lwd = 1)
text(x = sc_bp[c(4,5),1]*40 + 0.7,
     y = sc_bp[c(4,5),2]*40, 
     labels = rownames(sc_bp)[c(4,5)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(3),1]*40,
     y = sc_bp[c(3),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.2,
     y = sc_bp[c(1),2]*40 + 0.2, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40 + 0.5,
     y = sc_bp[c(2),2]*40 - 0.25, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
legend("bottomright", legend = c("Individual", "Gene"), 
       pch = c(21, 22), pt.bg = c(color_ind, color_gene), col = "black")
dev.off()

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

# --- custom RDA plot again, but with overlay of ge and twas genes
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1_de_twas.png", width = 850, height = 800)
plot(bmi_anc_rda,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlim = c(-3,3), 
     ylim = c(-3,3),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)"), 
     cex.lab = 1.5,
     cex.axis = 1.2
)
color_ind <- rgb(70/255, 130/255, 180/255, alpha = 0.5)  # steelblue with alpha
color_gene <- rgb(242/255, 189/255, 51/255, alpha = 0.5)
color_de <- "red"
color_twas <- "maroon"
points(sc_ind, 
       pch = 21,
       col = "white",
       bg = color_ind,
       cex = 1)
w_de <- which(rownames(sc_gene) %in% de_ensembl$ensembl_gene_id)
w_twas <- which(rownames(sc_gene) %in% genes_sig)
points(sc_gene[-c(w_de, w_twas),], 
       pch = 22,
       col = "white",
       bg = color_gene, 
       cex = 1)
points(sc_gene[c(w_twas),], 
       pch = 23,
       col = "white",
       bg = color_twas, 
       cex = 1.5)
points(sc_gene[c(w_de),], 
       pch = 24,
       col = "white",
       bg = color_de, 
       cex = 1.5)
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*40, sc_bp[,2]*40, # scale up the score values
       col = "red", 
       lwd = 1)
text(x = sc_bp[c(4,5),1]*40 + 0.7,
     y = sc_bp[c(4,5),2]*40, 
     labels = rownames(sc_bp)[c(4,5)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(3),1]*40,
     y = sc_bp[c(3),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.2,
     y = sc_bp[c(1),2]*40 + 0.2, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40 + 0.5,
     y = sc_bp[c(2),2]*40 - 0.25, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
legend("bottomright", legend = c("Individual", "Genes", "Known DEG", "Known TWAS"), 
       pch = c(21, 22, 23, 24), pt.bg = c(color_ind, color_gene, color_de, color_twas), col = "black")
dev.off()

# --- custom RDA triplot with de genes overlay, no twas
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1_de.png", width = 800, height = 800)
plot(bmi_anc_rda,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlim = c(-3,3), 
     ylim = c(-3,3),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)"), 
     cex.lab = 1.5,
     cex.axis = 1.2
)
color_ind <- rgb(70/255, 130/255, 180/255, alpha = 0.5)
color_gene <- rgb(242/255, 189/255, 51/255, alpha = 0.5)
color_de <- "red"
points(sc_ind, 
       pch = 21,
       col = "white",
       bg = color_ind,
       cex = 1)
w_de <- which(rownames(sc_gene) %in% de_ensembl$ensembl_gene_id)
w_twas <- which(rownames(sc_gene) %in% genes_sig)
points(sc_gene[-c(w_de),], 
       pch = 22,
       col = "white",
       bg = color_gene, 
       cex = 1)
points(sc_gene[c(w_de),], 
       pch = 24,
       col = "white",
       bg = color_de, 
       cex = 1.5)
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*40, sc_bp[,2]*40, # scale up the score values
       col = "red", 
       lwd = 1)
text(x = sc_bp[c(4,5),1]*40 + 0.7,
     y = sc_bp[c(4,5),2]*40, 
     labels = rownames(sc_bp)[c(4,5)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(3),1]*40,
     y = sc_bp[c(3),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.2,
     y = sc_bp[c(1),2]*40 + 0.2, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40 + 0.5,
     y = sc_bp[c(2),2]*40 - 0.25, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
legend("bottomright", legend = c("Individual", "Genes", "Known DEG"), 
       pch = c(21, 22, 23, 24), pt.bg = c(color_ind, color_gene, color_de), col = "black")
dev.off()


# --- loadings = where the new RDA axis lie with respect to the original axis 
gene_scores <- bmi_anc_rda$CCA$v
axis_loadings <- bmi_anc_rda$CCA$biplot

# --- can we quantify the gene's scores on each of the original axis? 
# gene_score_axis1 = gene_score_RDA1 * axis1_loading_RDA1 + gene_score_RDA2 * axis1_loading_RDA2
# (genes x RDA) %*% (RDA x axis)

axis_loadings_t <- t(axis_loadings)
gene_scores_axis <- gene_scores %*% axis_loadings_t

# --- order based on these gene axis score

gene_scores_axis_sorted_bmi <- gene_scores_axis[order(abs(gene_scores_axis[, "bmi"]), decreasing = TRUE), ] # sort based on absolute value of gene bmi score
genes <- rownames(gene_scores_axis_sorted_bmi)
gene_scores_axis_sorted_bmi_df <- as.data.frame(gene_scores_axis_sorted_bmi)
gene_scores_axis_sorted_bmi_df <- cbind(genes, gene_scores_axis_sorted_bmi_df)
rownames(gene_scores_axis_sorted_bmi_df) <- NULL
write.table(gene_scores_axis_sorted_bmi_df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/rda_bmi_sorted_genes.txt", quote = F, row.names = F, col.names = T, sep = '\t')

gene_scores_axis_sorted_CAUbmi_df <- gene_scores_axis_sorted_bmi_df[order(abs(gene_scores_axis_sorted_bmi_df[, "ancestryCAU:bmi"]), decreasing = TRUE), ] 
write.table(gene_scores_axis_sorted_CAUbmi_df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/rda_CAUbmi_sorted_genes.txt", quote = F, row.names = F, col.names = T, sep = '\t')

gene_scores_axis_sorted_HISbmi_df <- gene_scores_axis_sorted_bmi_df[order(abs(gene_scores_axis_sorted_bmi_df[, "ancestryHIS:bmi"]), decreasing = TRUE), ] 
write.table(gene_scores_axis_sorted_HISbmi_df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/rda_HISbmi_sorted_genes.txt", quote = F, row.names = F, col.names = T, sep = '\t')

# --- distribution of gene axis scores and plotting on original triplots
hist(gene_scores_axis_sorted_bmi_df$bmi) # gene-bmi scores normally distributed 
hist(gene_scores_axis_sorted_bmi_df$`ancestryCAU:bmi`) # gene-bmi*cau scores normally distributed 
hist(gene_scores_axis_sorted_bmi_df$`ancestryHIS:bmi`) # gene-bmi*his scores normally distributed 

top1_perc_bmi_genes <- gene_scores_axis_sorted_bmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_bmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmi_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #42.857% of the previously found DE genes are in the top 1% of BMI scored genes

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1_top1perc_bmi.png", width = 850, height = 800)
plot(bmi_anc_rda,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlim = c(-3,3), 
     ylim = c(-3,3),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)"), 
     cex.lab = 1.5,
     cex.axis = 1.2
)
color_ind <- rgb(70/255, 130/255, 180/255, alpha = 0.5)  # steelblue with alpha
color_gene <- rgb(242/255, 189/255, 51/255, alpha = 0.5)
color_bmi <- "maroon"
points(sc_ind, 
       pch = 21,
       col = "white",
       bg = color_ind,
       cex = 1)
w_bmi <- which(rownames(sc_gene) %in% top1_perc_bmi_genes$genes)
points(sc_gene[-c(w_de, w_twas),], 
       pch = 22,
       col = "white",
       bg = color_gene, 
       cex = 1)
points(sc_gene[c(w_bmi),], 
       pch = 23,
       col = "white",
       bg = color_bmi, 
       cex = 1.5)
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*40, sc_bp[,2]*40, # scale up the score values
       col = "red", 
       lwd = 1)
text(x = sc_bp[c(4,5),1]*40 + 0.7,
     y = sc_bp[c(4,5),2]*40, 
     labels = rownames(sc_bp)[c(4,5)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(3),1]*40,
     y = sc_bp[c(3),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.2,
     y = sc_bp[c(1),2]*40 + 0.2, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40 + 0.5,
     y = sc_bp[c(2),2]*40 - 0.25, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
legend("bottomright", legend = c("Individual", "Genes", "Top 1% BMI genes"), 
       pch = c(21, 22, 23), pt.bg = c(color_ind, color_gene, color_bmi), col = "black")
dev.off()

top1_perc_bmiHIS_genes <- gene_scores_axis_sorted_HISbmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_HISbmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmiHIS_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #16.7% are differentially expressed

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1_top1perc_bmiHIS.png", width = 850, height = 800)
plot(bmi_anc_rda,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlim = c(-3,3), 
     ylim = c(-3,3),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)"), 
     cex.lab = 1.5,
     cex.axis = 1.2
)
color_ind <- rgb(70/255, 130/255, 180/255, alpha = 0.5)  # steelblue with alpha
color_gene <- rgb(242/255, 189/255, 51/255, alpha = 0.5)
color_bmi <- "maroon"
points(sc_ind, 
       pch = 21,
       col = "white",
       bg = color_ind,
       cex = 1)
w_bmi <- which(rownames(sc_gene) %in% top1_perc_bmiHIS_genes$genes)
points(sc_gene[-c(w_de, w_twas),], 
       pch = 22,
       col = "white",
       bg = color_gene, 
       cex = 1)
points(sc_gene[c(w_bmi),], 
       pch = 23,
       col = "white",
       bg = color_bmi, 
       cex = 1.5)
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*40, sc_bp[,2]*40, # scale up the score values
       col = "red", 
       lwd = 1)
text(x = sc_bp[c(4,5),1]*40 + 0.7,
     y = sc_bp[c(4,5),2]*40, 
     labels = rownames(sc_bp)[c(4,5)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(3),1]*40,
     y = sc_bp[c(3),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.2,
     y = sc_bp[c(1),2]*40 + 0.2, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40 + 0.5,
     y = sc_bp[c(2),2]*40 - 0.25, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
legend("bottomright", legend = c("Individual", "Genes", "Top 1% BMI:HIS genes"), 
       pch = c(21, 22, 23), pt.bg = c(color_ind, color_gene, color_bmi), col = "black")
dev.off()

# --- cor and overlap of BMI scores and BMI:HIS scores

cor.test(gene_scores_axis_sorted_bmi_df$bmi, gene_scores_axis_sorted_bmi_df$`ancestryHIS:bmi`)$p.val
  # p < 1e-200, r = 0.5167016
ggplot(gene_scores_axis_sorted_bmi_df, aes(x = bmi, y = `ancestryHIS:bmi`)) +
  geom_point() +  # Add points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  geom_smooth(method = "lm", color = "blue", se = FALSE) + 
  labs(x = "Gene Scores on BMI axis", y = "Gene Scores on BMI:HIS axis") + 
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/bmi_bmiHIS_cor.png", width = 5, height = 5)

top1_perc_bmiHIS_genes
top1_perc_bmi_genes

# --- BMI vs ancestry
boxplot(bmi ~ ancestry, type = "n") # seems to be some difference in bmi across ancestries
bmi_anc <- lm(bmi ~ ancestry)
bmi_anc_anova <- anova(bmi_anc)
bmi_anc_lm <- summary(bmi_anc)
df_plot <- as.data.frame(cbind(ancestry, as.numeric(bmi)))
ggplot(df_plot, aes(x = ancestry, y = bmi, fill = ancestry)) + 
  geom_boxplot() +
  labs(x = "Ancestry", y = "BMI") +
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/BMI_anc.png", width = 6, height = 6)

