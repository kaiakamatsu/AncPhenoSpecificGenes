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


# --- BMI vs ancestry, preliminary analysis 
ancestry <- df$race
bmi <- df$`bmi_(kg/m^2)`
boxplot(bmi ~ ancestry, type = "n") # seems to be some difference in bmi across ancestries
bmi_anc <- lm(bmi ~ ancestry)
bmi_anc_anova <- anova(bmi_anc)
TukeyHSD(aov(bmi ~ ancestry))
bmi_anc_lm <- summary(bmi_anc)
df_plot <- as.data.frame(cbind(ancestry, as.numeric(bmi)))
ggplot(df_plot, aes(x = ancestry, y = bmi, fill = ancestry)) + 
  geom_boxplot() +
  labs(x = "Ancestry", y = "BMI") +
  theme_bw() +
  theme_classic(base_size = 16)
ggsave(filename = "C:/Users/kaiak/OneDrive/Documents/BGGN273/final_project/AncPhenoSpecificGenes/BMI_anc.png", width = 6, height = 6)


