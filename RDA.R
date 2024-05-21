install.packages('vegan')
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(fBasics)
library(MASS)

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

# --- plot RDA results 
png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_vegan_scale1.png")
ordiplot(bmi_anc_rda, scaling = 1, type = "points")
dev.off()
# interpretation is difficult, constrained variance is very small so we are only looking at 1.3% of the total variance in this plot

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_vegan_scale2.png")
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

png(filename="C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/RDA_custom_scale1.png")
plot(bmi_anc_rda,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlim = c(-3,3), 
     ylim = c(-3,3),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
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
text(x = sc_bp[c(4,5),1]*40 + 1,
     y = sc_bp[c(4,5),2]*40, 
     labels = rownames(sc_bp)[c(4,5)], 
     col = "red", 
     cex = 1, 
     font = 2)
text(x = sc_bp[c(3),1]*40,
     y = sc_bp[c(3),2]*40 - 0.2, 
     labels = rownames(sc_bp)[c(3)], 
     col = "red", 
     cex = 1, 
     font = 2)
text(x = sc_bp[c(1),1]*40 + 0.2,
     y = sc_bp[c(1),2]*40 + 0.2, 
     labels = rownames(sc_bp)[c(1)], 
     col = "red", 
     cex = 1, 
     font = 2)
text(x = sc_bp[c(2),1]*40 + 0.5,
     y = sc_bp[c(2),2]*40 - 0.25, 
     labels = rownames(sc_bp)[c(2)], 
     col = "red", 
     cex = 1, 
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
pc_scores <- pca$x
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = "maroon")+ 
  theme_bw() +
  theme_classic(base_size = 14)



# --- BMI vs ancestry
boxplot(bmi ~ ancestry, type = "n") # seems to be some difference in bmi across ancestries
bmi_anc <- lm(bmi ~ ancestry)
bmi_anc_anova <- anova(bmi_anc)
bmi_anc_lm <- summary(bmi_anc)


