library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(fBasics)
library(MASS)
library(biomaRt)

### bmi_anc_rda is the rda object from vegan

# --- plot RDA results 
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_vegan_scale1.png", width = 800, height = 800)
ordiplot(bmi_anc_rda, scaling = 1, type = "points")
dev.off()

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

# --- TWAS BMI genes 
twas <- fread("C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/twas.dat")
unique(twas$PANEL) # a tone of tissue panels 
# these expression values are in Monocytes, so work with blood related panels? 
twas_blood <- filter(twas, PANEL == "GTEx.Whole_Blood")
w <- which(twas_blood$TWAS.P < (0.05/nrow(twas_blood))) # bonferroni
twas_blood_significant <- twas_blood[w,]
genes_sig <- sapply(twas_blood_significant$FILE, function(x) strsplit(x, split = "[.]")[[1]][5])


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



