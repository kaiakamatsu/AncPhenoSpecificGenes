library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(fBasics)
library(MASS)
library(biomaRt)


# --- read in data again 
ge_residuals <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/processed_ge_after_covar.txt")
df <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/transformed_phenos.txt")
covar <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/covars.txt")

# --- RDA model
ancestry <- df$race
bmi <- df$`bmi_(kg/m^2)`
bmi_anc_rda <- rda(formula = ge_residuals ~ ancestry + bmi)

# --- hypothesis test
null_rda <- rda(formula = ge_residuals ~ 1)
test <- anova.cca(null_rda, bmi_anc_rda)
#Model 1: ge_residuals ~ 1
#Model 2: ge_residuals ~ ancestry + bmi
#ResDf ResChiSquare Df ChiSquare      F Pr(>F)    
#1  1110        22475                               
#2  1107        22209  3    265.81 4.4164  0.001 ***

# --- custom RDA plot 
perc <- round(100*(summary(bmi_anc_rda)$cont$importance[2, 1:2]), 2)
sc_ind <- scores(bmi_anc_rda, display="sites", choices=c(1,2), scaling=1)
sc_gene <- scores(bmi_anc_rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(bmi_anc_rda, display="bp", choices=c(1, 2), scaling=1)

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
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/BMI_RDA_custom_scale1_de.png", width = 800, height = 800)
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
text(x = sc_bp[c(3),1]*40 + 0.2,
     y = sc_bp[c(3),2]*40, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.3,
     y = sc_bp[c(1),2]*40 + 0.3, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40,
     y = sc_bp[c(2),2]*40 - 0.2, 
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
axis_loadings_t <- t(axis_loadings)*eigen_values # axis loading scaled by eigen values
gene_scores_axis <- gene_scores %*% axis_loadings_t

# --- order based on these gene axis score
gene_scores_axis_sorted_bmi <- gene_scores_axis[order(abs(gene_scores_axis[, "bmi"]), decreasing = TRUE), ] # sort based on absolute value of gene bmi score
genes <- rownames(gene_scores_axis_sorted_bmi)
gene_scores_axis_sorted_bmi_df <- as.data.frame(gene_scores_axis_sorted_bmi)
gene_scores_axis_sorted_bmi_df <- cbind(genes, gene_scores_axis_sorted_bmi_df)
rownames(gene_scores_axis_sorted_bmi_df) <- NULL
write.table(gene_scores_axis_sorted_bmi_df, file = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/rda_bmi_sorted_genes.txt", quote = F, row.names = F, col.names = T, sep = '\t')

# --- distribution of gene axis scores and plotting on original triplots
hist(gene_scores_axis_sorted_bmi_df$bmi) # gene-bmi scores normally distributed 
top1_perc_bmi_genes <- gene_scores_axis_sorted_bmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_bmi_df)/100)) ),]
#top10_perc_bmi_genes <- gene_scores_axis_sorted_bmi_df[c(1:as.integer((nrow(gene_scores_axis_sorted_bmi_df)/10)) ),]

w_overlap_de <- which(top1_perc_bmi_genes$genes %in% de_ensembl$ensembl_gene_id)
(length(w_overlap_de)/nrow(de_ensembl))*100 #19% of the previously found DE genes are in the top 1% of BMI scored genes
  # 8/42 differentially expressed genes 


# --- visualize top 1 percent, BMI overlap
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/BMI_RDA_custom_scale1_top1perc_bmi.png", width = 850, height = 800)
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
color_bmi <- "darkgreen"
color_de <- "darkblue"
color_both <- "maroon"
points(sc_ind, 
       pch = 21,
       col = "white",
       bg = color_ind,
       cex = 1)
w_both <- which(rownames(sc_gene) %in% de_ensembl$ensembl_gene_id & rownames(sc_gene) %in% top1_perc_bmi_genes$genes)
w_de <- setdiff(which(rownames(sc_gene) %in% de_ensembl$ensembl_gene_id), w_both)
w_bmi <- setdiff(which(rownames(sc_gene) %in% top1_perc_bmi_genes$genes), w_both)
points(sc_gene[-c(w_de, w_bmi, w_both),], 
       pch = 22,
       col = "white",
       bg = color_gene, 
       cex = 1)
points(sc_gene[c(w_bmi),], 
       pch = 22,
       col = "white",
       bg = color_bmi, 
       cex = 2)
points(sc_gene[c(w_de),], 
       pch = 22,
       col = "white",
       bg = color_de, 
       cex = 2)
points(sc_gene[c(w_both),], 
       pch = 23,
       col = "white",
       bg = color_both, 
       cex = 3)
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*40, sc_bp[,2]*40, # scale up the score values
       col = "red", 
       lwd = 1)
text(x = sc_bp[c(3),1]*40 + 0.2,
     y = sc_bp[c(3),2]*40, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.3,
     y = sc_bp[c(1),2]*40 + 0.3, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
text(x = sc_bp[c(2),1]*40,
     y = sc_bp[c(2),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1.2, 
     font = 2)
legend("bottomright", legend = c("Individual", "Genes", "DE", "Top 1% BMI genes", "Both DE and BMI"), 
       pch = c(21, 22, 22, 22, 23), pt.bg = c(color_ind, color_gene, color_de, color_bmi, color_both), col = "black")
dev.off()


# --- top 1 percent of bmi genes
w_cols_bmi <- which(colnames(ge_residuals) %in% top1_perc_bmi_genes$genes)
ge_top1_perc_bmi <- as.data.frame(ge_residuals)[,w_cols_bmi]

bmi_effects <- function(ge_values, bmi_values){
  model <- lm(ge_values ~ bmi_values)
  model_summary <- summary(model)
  data <- model_summary$coefficients[2,c(1,2,4)]
  return(data)
}

genes <- c()
B <- c()
SE <- c()
p <- c()
for(i in 1:ncol(ge_top1_perc_bmi)){
  genes <- append(genes, colnames(ge_top1_perc_bmi)[i])
  datas <- bmi_effects(ge_top1_perc_bmi[,i], bmi)
  B <- append(B, datas[1])
  SE <- append(SE, datas[2])
  p <- append(p, datas[3])
}
nrow(summary_stats) # 224 total tests
summary_stats <- data.frame(genes, B, SE, p)
summary_stats_sig <- filter(summary_stats, p < (0.05/nrow(gene_scores_axis_sorted_bmi_df)))
summary_stats_sig_lenient <- filter(summary_stats, p < (0.05/ (nrow(summary_stats)*3) ))
nrow(summary_stats_sig) # 144 hits, bonferroni 
nrow(summary_stats_sig_lenient) # 193 hits, bonferroni 

# --- plot the strongest association 
strongest_association = summary_stats[which(summary_stats$p == min(summary_stats$p) ),]
top_gene_expression <- ge_top1_perc_bmi[,which(colnames(ge_top1_perc_bmi) == strongest_association$genes)]
lm_best <- lm(top_gene_expression ~ bmi + ancestry)
step(lm_best)
sum <- summary(lm_best) # intercepts are significant
df_plot_top <- data.frame(bmi, ancestry, top_gene_expression)

cc <- c("maroon", "darkblue", "darkgreen")
ggplot(df_plot_top, aes(x = bmi, y = top_gene_expression, color = ancestry)) +
  geom_point() + 
  scale_color_manual(values = cc) +  
  labs(x = "BMI", y = "RRAGD expression") + 
  theme_bw() +
  theme_classic(base_size = 16) +
  geom_abline(slope = 0.2284, intercept = 0.4092, color = "maroon", linewidth = 1.3) +  
  geom_abline(slope = 0.2284, intercept = (0.4092 - 0.5131), color = "darkblue", linewidth = 1.3) + 
  geom_abline(slope = 0.2284, intercept = (0.4092 - 0.4894), color = "darkgreen", linewidth = 1.3)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/figures/RRAGDvsBMI.png", width = 6, height = 5)


