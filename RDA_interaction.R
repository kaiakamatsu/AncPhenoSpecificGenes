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
eigen_values <- bmi_anc_rda$CCA$eig

# --- can we quantify the gene's scores on each of the original axis?
# gene_score_axis1 = gene_score_RDA1 * axis1_loading_RDA1 + gene_score_RDA2 * axis1_loading_RDA2
# (genes x RDA) %*% (RDA x axis)

axis_loadings_t <- t(axis_loadings)*eigen_values
gene_scores_axis <- gene_scores %*% axis_loadings_t

# --- order based on these gene axis score

gene_scores_axis_sorted_CAUbmi <- gene_scores_axis[order(abs(gene_scores_axis[, "ancestryCAU:bmi"]), decreasing = TRUE), ]
genes <- rownames(gene_scores_axis_sorted_CAUbmi)
gene_scores_axis_sorted_CAUbmi_df <- as.data.frame(gene_scores_axis_sorted_CAUbmi)
gene_scores_axis_sorted_CAUbmi_df <- cbind(genes, gene_scores_axis_sorted_CAUbmi_df)
rownames(gene_scores_axis_sorted_CAUbmi_df) <- NULL
write.table(gene_scores_axis_sorted_CAUbmi_df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/rda_CAUbmi_sorted_genes.txt", quote = F, row.names = F, col.names = T, sep = '\t')

gene_scores_axis_sorted_HISbmi_df <- gene_scores_axis_sorted_CAUbmi_df[order(abs(gene_scores_axis_sorted_CAUbmi_df[, "ancestryHIS:bmi"]), decreasing = TRUE), ] 
write.table(gene_scores_axis_sorted_HISbmi_df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/rda_HISbmi_sorted_genes.txt", quote = F, row.names = F, col.names = T, sep = '\t')

gene_scores_axis_sorted_bmi_df <- gene_scores_axis_sorted_CAUbmi_df[order(abs(gene_scores_axis_sorted_CAUbmi_df[, "bmi"]), decreasing = TRUE), ] 


# --- top 1 percent bmi:ANC genes
top1_perc_bmiHIS_genes <- gene_scores_axis_sorted_HISbmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_HISbmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmiHIS_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #30% are differentially expressed

top1_perc_bmiCAU_genes <- gene_scores_axis_sorted_CAUbmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_CAUbmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmiCAU_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #21%

top1_perc_bmi_genes <- gene_scores_axis_sorted_bmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_bmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmi_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #19%

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
points(sc_gene[-c(w_bmi),], 
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

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1_top1perc_bmiCAU.png", width = 850, height = 800)
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
w_bmi <- which(rownames(sc_gene) %in% top1_perc_bmiCAU_genes$genes)
points(sc_gene[-c(w_bmi),], 
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
legend("bottomright", legend = c("Individual", "Genes", "Top 1% BMI:CAU genes"), 
       pch = c(21, 22, 23), pt.bg = c(color_ind, color_gene, color_bmi), col = "black")
dev.off()

# --- cor and overlap of BMI scores and BMI:HIS scores
cor.test(gene_scores_axis_sorted_bmi_df$bmi, gene_scores_axis_sorted_bmi_df$`ancestryHIS:bmi`)$p.val
  # p < 1e-200, r = 0.8441867
ggplot(gene_scores_axis_sorted_bmi_df, aes(x = bmi, y = `ancestryHIS:bmi`)) +
  geom_point() +  # Add points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  geom_smooth(method = "lm", color = "blue", se = FALSE) + 
  labs(x = "Gene Scores on BMI axis", y = "Gene Scores on BMI:HIS axis") + 
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/bmi_bmiHIS_cor.png", width = 5, height = 5)

cor.test(gene_scores_axis_sorted_bmi_df$bmi, gene_scores_axis_sorted_bmi_df$`ancestryCAU:bmi`)$p.val
# p < 1e-200, r = 0.9390707
ggplot(gene_scores_axis_sorted_bmi_df, aes(x = bmi, y = `ancestryCAU:bmi`)) +
  geom_point() +  # Add points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  geom_smooth(method = "lm", color = "blue", se = FALSE) + 
  labs(x = "Gene Scores on BMI axis", y = "Gene Scores on BMI:CAU axis") + 
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/bmi_bmiCAU_cor.png", width = 5, height = 5)


intersection <- length(which(top1_perc_bmi_genes$genes %in% top1_perc_bmiHIS_genes$genes)) #126 genes in common
specific <- nrow(top1_perc_bmi_genes) - intersection #98 unique to each group

# --- top 1 percent of bmi:his genes 
w_cols_bmi_his <- which(colnames(ge_residuals) %in% top1_perc_bmiHIS_genes$genes)
ge_top1_perc_bmi_his <- as.data.frame(ge_residuals)[,w_cols_bmi_his]

bmi_anc_interaction_effects <- function(ge_values, bmi_values, ancestry_groups, col){
  model <- lm(ge_values ~ bmi_values*ancestry_groups)
  model_summary <- summary(model)
  data <- model_summary$coefficients[col,c(4)]
  return(data)
}

genes <- c()
p_HIS <- c()
for(i in 1:ncol(ge_top1_perc_bmi_his)){
  genes <- append(genes, colnames(ge_top1_perc_bmi_his)[i])
  datas <- bmi_anc_interaction_effects(ge_top1_perc_bmi_his[,i], bmi, ancestry, 6)
  p_HIS <- append(p_HIS, datas)
}
summary_stats_bmianc <- data.frame(genes, p_HIS)
summary_stats_bmiHIS_sig <- filter(summary_stats_bmianc, p_HIS < (0.05/ (nrow(summary_stats_bmianc)*5) ))
nrow(summary_stats_bmiHIS_sig)


# --- top 1 percent of bmi:cau genes 
w_cols_bmi_cau <- which(colnames(ge_residuals) %in% top1_perc_bmiCAU_genes$genes)
ge_top1_perc_bmi_cau <- as.data.frame(ge_residuals)[,w_cols_bmi_cau]

genes <- c()
p_CAU <- c()
for(i in 1:ncol(ge_top1_perc_bmi_his)){
  genes <- append(genes, colnames(ge_top1_perc_bmi_cau)[i])
  datas <- bmi_anc_interaction_effects(ge_top1_perc_bmi_cau[,i], bmi, ancestry, 5)
  p_CAU <- append(p_CAU, datas)
}
summary_stats_bmicau <- data.frame(genes, p_CAU)
summary_stats_bmiCAU_sig <- filter(summary_stats_bmicau, p_CAU < (0.05/ (nrow(summary_stats_bmianc)*5) ))
nrow(summary_stats_bmiCAU_sig) 






### --- NOT WEIGHTED BY EIGENVALUES 

# --- loadings = where the new RDA axis lie with respect to the original axis 
gene_scores <- bmi_anc_rda$CCA$v
axis_loadings <- bmi_anc_rda$CCA$biplot

# --- can we quantify the gene's scores on each of the original axis?
# gene_score_axis1 = gene_score_RDA1 * axis1_loading_RDA1 + gene_score_RDA2 * axis1_loading_RDA2
# (genes x RDA) %*% (RDA x axis)
axis_loadings_t <- t(axis_loadings)
gene_scores_axis <- gene_scores %*% axis_loadings_t

# --- order based on these gene axis score
gene_scores_axis_sorted_CAUbmi <- gene_scores_axis[order(abs(gene_scores_axis[, "ancestryCAU:bmi"]), decreasing = TRUE), ]
genes <- rownames(gene_scores_axis_sorted_CAUbmi)
gene_scores_axis_sorted_CAUbmi_df <- as.data.frame(gene_scores_axis_sorted_CAUbmi)
gene_scores_axis_sorted_CAUbmi_df <- cbind(genes, gene_scores_axis_sorted_CAUbmi_df)
rownames(gene_scores_axis_sorted_CAUbmi_df) <- NULL

gene_scores_axis_sorted_HISbmi_df <- gene_scores_axis_sorted_CAUbmi_df[order(abs(gene_scores_axis_sorted_CAUbmi_df[, "ancestryHIS:bmi"]), decreasing = TRUE), ] 

gene_scores_axis_sorted_bmi_df <- gene_scores_axis_sorted_CAUbmi_df[order(abs(gene_scores_axis_sorted_CAUbmi_df[, "bmi"]), decreasing = TRUE), ] 

# --- overlap to prev found de genes
top1_perc_bmiHIS_genes <- gene_scores_axis_sorted_HISbmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_HISbmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmiHIS_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #17% are differentially expressed

top1_perc_bmiCAU_genes <- gene_scores_axis_sorted_CAUbmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_CAUbmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmiCAU_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #17%

top1_perc_bmi_genes <- gene_scores_axis_sorted_bmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_bmi_df)/100)) ),]
w_overlap_de <- which(top1_perc_bmi_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #43%

# --- top 1 percent of bmi:his genes 
bmi_anc_interaction_effects <- function(ge_values, bmi_values, ancestry_groups, col){
  model <- lm(ge_values ~ bmi_values*ancestry_groups)
  model_summary <- summary(model)
  data <- model_summary$coefficients[col,c(4)]
  return(data)
}

w_cols_bmi_his <- which(colnames(ge_residuals) %in% top1_perc_bmiHIS_genes$genes)
ge_top1_perc_bmi_his <- as.data.frame(ge_residuals)[,w_cols_bmi_his]

genes <- c()
p_HIS <- c()
for(i in 1:ncol(ge_top1_perc_bmi_his)){
  genes <- append(genes, colnames(ge_top1_perc_bmi_his)[i])
  datas <- bmi_anc_interaction_effects(ge_top1_perc_bmi_his[,i], bmi, ancestry, 6)
  p_HIS <- append(p_HIS, datas)
}
summary_stats_bmianc <- data.frame(genes, p_HIS)
summary_stats_bmiHIS_sig <- filter(summary_stats_bmianc, p_HIS < (0.05/ (nrow(summary_stats_bmianc)*5) ))
nrow(summary_stats_bmiHIS_sig)

# --- top 1 percent of bmi:cau genes 
w_cols_bmi_cau <- which(colnames(ge_residuals) %in% top1_perc_bmiCAU_genes$genes)
ge_top1_perc_bmi_cau <- as.data.frame(ge_residuals)[,w_cols_bmi_cau]

genes <- c()
p_CAU <- c()
for(i in 1:ncol(ge_top1_perc_bmi_his)){
  genes <- append(genes, colnames(ge_top1_perc_bmi_cau)[i])
  datas <- bmi_anc_interaction_effects(ge_top1_perc_bmi_cau[,i], bmi, ancestry, 5)
  p_CAU <- append(p_CAU, datas)
}
summary_stats_bmicau <- data.frame(genes, p_CAU)
summary_stats_bmiCAU_sig <- filter(summary_stats_bmicau, p_CAU < (0.05/ (nrow(summary_stats_bmianc)*5) ))
nrow(summary_stats_bmiCAU_sig) 

# --- plot
strongest_association = summary_stats_bmicau[which(summary_stats_bmicau$p == min(summary_stats_bmicau$p) ),]
top_gene_expression <- ge_top1_perc_bmi_cau[,which(colnames(ge_top1_perc_bmi_cau) == strongest_association$genes)]
lm_best <- lm(top_gene_expression ~ bmi*ancestry)
step(lm_best)
sum <- summary(lm_best) # intercepts not significant, slopes diff, split the datasets
df_plot_top <- data.frame(bmi, ancestry, top_gene_expression)

df_plot_top_aa <- df_plot_top[which(df_plot_top$ancestry == "AA"),]
lm_aa <- lm(df_plot_top_aa$top_gene_expression ~ df_plot_top_aa$bmi)
step(lm_aa)
summary(lm_aa) 
df_plot_top_his <- df_plot_top[which(df_plot_top$ancestry == "HIS"),]
lm_his <- lm(df_plot_top_his$top_gene_expression ~ df_plot_top_his$bmi)
step(lm_his)
summary(lm_his)
df_plot_top_cau <- df_plot_top[which(df_plot_top$ancestry == "CAU"),]
lm_cau <- lm(df_plot_top_cau$top_gene_expression ~ df_plot_top_cau$bmi)
step(lm_cau)
summary(lm_cau)


cc <- c("maroon", "darkblue", "darkgreen")
ggplot(df_plot_top, aes(x = bmi, y = top_gene_expression, color = ancestry)) +
  geom_point() + 
  scale_color_manual(values = cc) +  
  labs(x = "BMI", y = "TBCK expression") + 
  theme_bw() +
  theme_classic(base_size = 16) +
  geom_abline(slope = -0.210294, intercept = 0.040296, color = "maroon", linewidth = 1.3) +  
  geom_abline(slope = (-0.210294 + 0.397412), intercept = 0.040296, color = "darkblue", linewidth = 1.3) + 
  geom_abline(slope = (-0.210294 + 0.298217), intercept = 0.040296, color = "darkgreen", linewidth = 1.3)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/TBCKvsBMI.png", width = 6, height = 5)

cc <- c("maroon")
ggplot(df_plot_top_aa, aes(x = bmi, y = top_gene_expression, color = ancestry)) +
  geom_point() + 
  scale_color_manual(values = cc) +  
  labs(x = "BMI", y = "TBCK expression") + 
  theme_bw() +
  theme_classic(base_size = 16) +
  geom_abline(slope = -0.21029, intercept = 0.04030, color = "maroon", linewidth = 1.3)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/TBCKvsBMI_AA.png", width = 6, height = 5)

cc <- c("darkblue")
ggplot(df_plot_top_cau, aes(x = bmi, y = top_gene_expression, color = ancestry)) +
  geom_point() + 
  scale_color_manual(values = cc) +  
  labs(x = "BMI", y = "TBCK expression") + 
  theme_bw() +
  theme_classic(base_size = 16) +
  geom_abline(slope = 0.18712, intercept = 0.04888, color = "darkblue", linewidth = 1.3)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/TBCKvsBMI_CAU.png", width = 6, height = 5)

cc <- c("darkgreen")
ggplot(df_plot_top_his, aes(x = bmi, y = top_gene_expression, color = ancestry)) +
  geom_point() + 
  scale_color_manual(values = cc) +  
  labs(x = "BMI", y = "TBCK expression") + 
  theme_bw() +
  theme_classic(base_size = 16) +
  geom_abline(slope = 0.08792, intercept = -0.01700, color = "darkgreen", linewidth = 1.3)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/TBCKvsBMI_HIS.png", width = 6, height = 5)


