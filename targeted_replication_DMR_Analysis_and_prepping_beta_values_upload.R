###### DMR replication Analysis


########################################## set working directory and load the needed packages

## set workig directory
setwd('/Users/laniemullins/Desktop/MSPH_Thesis/Grady_Trauma_Project/')

## clear some space
rm(list=ls())

## load 
library(tibble)
if(!require("BiocManager", quietly = TRUE))
  install.packages(BiocManager)
BiocManager::install(version = "3.18")
BiocManager::install("mCSEA")
library(mCSEA)
library(xlsx)
library(ggpubr)
BiocManager::install(("clusterProfiler"))
library(clusterProfiler)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(data.table)
install.packages("dplyr")
library(dplyr)
library(biomaRt)
install.packages("org.Hs.eg.db")
library(org.Hs.eg.db)
install.packages("Hmisc")
library(Hmisc)
install.packages("BiocFileCache")

## set seed
set.seed(42)


########################## defining needed functions to calculate variance and
## select the probes with the top 5% var

## we will select the top cpgs based on variance
## First we need to calculated the variance for each row (cpg)

get_variance <- function(beta){
  # calculate variance for each row 
  beta$var <- apply(beta, 1, var)
  
  # order the data frame based on varaiences in descending order
  x <- beta[order(beta$var, decreasing = TRUE), ]
  
  # return ordered data frame
  return(x)
}

## we want top 5% probes with high variance
## select the top cpgs based on variance
get_top <- function(beta, per){
  
  # calculate the number of top rows based on the specificed percentage
  top <- (nrow(beta) * per)/100
  
  # select the top rows based on varience
  x <- beta[1:top, ]
  return(x)
}

############################ prep methylation values

## r prepare/format data never vs current
## merge files of interest into 3 pheno and beta files 
#### pheno_1/beta_1 is never and current, pheno_2/beta_2 is never, remitted, pheno_3/beta_3 is remitted, current

########################### running for never vs current
beta_vals_1 <- fread('GTP_Never_Current_ASM_Beta.csv')
head(beta_vals_1)
beta_vals_1 <- column_to_rownames(beta_vals_1, var = 'V1')
head(beta_vals_1)
#column_to_delete <- 1
#beta_vals_1 <- beta_vals_1[, -column_to_delete, drop = FALSE]

## convert beta values to a data frame to avoid the automatic conversion of character vectors to factors
beta_vals_1 <- as.data.frame(beta_vals_1, stringsAsFactors = FALSE)
dim(beta_vals_1)

# make sure all of the data is numeric, it not remove the column
# for some reason, there is an inserted column that is the probe names; remove it
data_types <- sapply(beta_vals_1, class)
print(data_types)
## no need to delete, all columns are numeric
#column_to_delete <- "V34"
#beta_vals_1 <- beta_vals_1[, !(colnames(beta_vals_1) %in% column_to_delete), drop = FALSE]

# Call functions, using only the top 5% probes
beta_1_wd_var <- get_variance(beta_vals_1)
# insure variences are not NA
head(beta_1_wd_var)

beta_1_top <- get_top(beta = beta_1_wd_var, per = 5)


dim(beta_1_top)
## 40976 beta values and 60 samples, as expected

## remove last column -- the variance from the DF
delete_var_column <- 61
beta_1_top <- beta_1_top[, -delete_var_column, drop = FALSE]
head(beta_1_top)

write.csv(beta_1_top, "top_beta_GTP_Never_Current_ASM.csv", quote = FALSE)

##prepare/format data never vs remitted

beta_vals_2 <- fread('GTP_Never_Remitted_ASM_Beta.csv')
head(beta_vals_2)
beta_vals_2 <- column_to_rownames(beta_vals_2, var = 'V1')
head(beta_vals_2)


## convert beta values to a data frame to avoid the automatic conversion of character vectors to factors
beta_vals_2 <- as.data.frame(beta_vals_2, stringsAsFactors = FALSE)
dim(beta_vals_2)

# make sure all of the data is numeric, it not remove the column
# for some reason, there is an inserted column that is the probe names; remove it
data_types <- sapply(beta_vals_2, class)
print(data_types)
## no need to delete, all columns are numeric


# Call functions, using only the top 5% probes
beta_2_wd_var <- get_variance(beta_vals_2)
# insure variences are not NA
head(beta_2_wd_var)

beta_2_top <- get_top(beta = beta_2_wd_var, per = 5)


dim(beta_2_top)
## 40976 beta values and 60 samples, as expected

## remove last column -- the variance from the DF
delete_var_column <- 83
beta_2_top <- beta_2_top[, -delete_var_column, drop = FALSE]
head(beta_1_top)

write.csv(beta_1_top, "top_beta_GTP_Never_Remitted_ASM.csv", quote = FALSE)


########################### running for never vs remitted
beta_vals_3 <- fread('GTP_Current_Remitted_ASM_Beta.csv')
head(beta_vals_3)
beta_vals_3 <- column_to_rownames(beta_vals_3, var = 'V1')
head(beta_vals_3)


## convert beta values to a data frame to avoid the automatic conversion of character vectors to factors
beta_vals_3 <- as.data.frame(beta_vals_3, stringsAsFactors = FALSE)
dim(beta_vals_3)

# make sure all of the data is numeric, it not remove the column
# for some reason, there is an inserted column that is the probe names; remove it
data_types <- sapply(beta_vals_3, class)
print(data_types)
## no need to delete, all columns are numeric


# Call functions, using only the top 5% probes
beta_3_wd_var <- get_variance(beta_vals_3)
# insure variences are not NA
head(beta_3_wd_var)

beta_3_top <- get_top(beta = beta_3_wd_var, per = 5)


dim(beta_3_top)
## 40976 beta values and 60 samples, as expected

## remove last column -- the variance from the DF
delete_var_column <- 72
beta_3_top <- beta_3_top[, -delete_var_column, drop = FALSE]
head(beta_3_top)

write.csv(beta_1_top, "top_beta_GTP_Current_Remitted_ASM.csv", quote = FALSE)


################################################ define functions for DMR analysis

## pull out required columns
get_required_cols <- function(pheno, cols){
  
  ## extract columns names using from the df
  x <- pheno[, cols]
  return(x)
}

## Run the test
## looks for DMRs in promoter regions
run_mCSEA <- function(b_rank, beta, pheno, annotation){
  x <- mCSEATest(b_rank, beta, pheno,
                 regionsTypes = "promoter", 
                 customAnnotation = annotation, 
                 platform = "EPIC",
                 nproc = 10, minCpGs = 5)
  return(x)
}



## select those with padj < 0.05
get_significant <- function(results){
  x <- results$promoters[results$promoters$padj < 0.05, ]
  return(x)
}

## function to select all which we want to do if we want 
## to have a background DMRs for plotting 

get_all_DMRs <- function(results){
  x <- results$promoters
  return(x)
} 

## function to plot the significant genes
## the output of mCSEA
## it will save pdfs
plot_genes <- function(path, results, genelis){
  setwd(path)
  for(i in 1:length(genelis)){
    g <- genelis[i]
    message('Saving plot for gene ', g)
    output <- tryCatch({
      mCSEAPlot(results, regionType = "promoters",
                dmrName = g, leadingEdge = TRUE,
                CGI = TRUE,
                transcriptAnnotation = "symbol", makePDF = TRUE)
    }, error = function(err){
      message("Gene not ploted ", g)
    })
  }
}


## turn cpgs from mCSEA into a vector
## function to split the cpgs of one gene

clean_cpgs <- function(cpgs){
  x <- unlist(strsplit(cpgs, ","))
  x <- gsub(" ", "", x, fixed = TRUE) # remove any spaces
  return(x)
}


## function to plot cpgs
plot_cpgs <- function(b_vals, g_name=NULL, s_name ){
  pl <- list()
  cpgs <- colnames(b_vals)[which(grepl('^cg', colnames(b_vals)))]
  for(i in 1:length(cpgs)){
    p <- ggboxplot(b_vals, x = "Exposure", y = cpgs[i],
                   color = "Exposure", palette = c("#0072B5FF", "#BC3C29FF"),# "nejm",
                   size = 1,
                   add = "dotplot"
    )
    ##  Add p-value
    p <- p + stat_compare_means(method ="t.test", label = "p.format") +
      font("title", size = 20)+ #face = "bold"
      font("xlab", size = 16)+
      font("ylab", size = 16)+
      font("xy.text", size = 13)+
      font("caption", size =15, face= "bold")+
      theme(legend.position = "none",
            axis.title.x =element_blank(),
            axis.title.y =element_blank())
    
    pl[[length(pl) + 1]] <- p
  }
  names(pl) <- cpgs
  return(pl)
}

##
get_sumry <- function(beta_df){
  print(dim(beta_df))
  print("Min: ")
  print(min(beta_df, na.rm = T))
  print("Max: ")
  print(max(beta_df, na.rm = T))
  col_mean <- colMeans(beta_df, na.rm = T)
  mean_ch <- sum(col_mean)/length(col_mean)
  print("mean: ")
  print(mean_ch)
  return(mean_ch)
}


############################################# freeing space

rm(beta_vals_1)
rm(beta_1_wd_var)
rm(beta_vals_2)
rm(beta_2_wd_var)
rm(beta_vals_3)
rm(beta_3_wd_var)



############################################ prepare pheno files

## pheno for never vs current}

## load phenotype data

#GTP_Never_Current_ASM_Beta.csv
#GTP_Never_Current_ASM_Pheno.csv

pheno_1 <- read.csv('GTP_Never_Current_ASM_Pheno.csv', 
                    header = TRUE, 
                    sep = ",")

pheno_1 <- column_to_rownames(pheno_1, var = 'SampleID')

column_to_delete <- 1
pheno_1 <- pheno_1[, -column_to_delete, drop = FALSE]
head(pheno_1)

# check that cell types sum to 1
## use the dplyr:: because R was confused about pathways
cell_sum <- pheno_1 %>% 
  dplyr::select(c("cd8t", "cd4t", "nk", "bcell","mono", "neu")) %>% 
  rowSums()
## all sums == 1

## Get only the information we need
#Binary MDD state, sex, cell type proportion
cols_of_interest <- c("explainatory", "cd8t", "cd4t", "nk", "bcell","mono", "neu", "Sex")
pheno_1_C <- pheno_1[cols_of_interest]

## prep pheno for never vs remitted}
## load phenotype data


pheno_2 <- read.csv('GTP_Never_Remitted_ASM_Pheno.csv', 
                    header = TRUE, 
                    sep = ",")

pheno_2 <- column_to_rownames(pheno_2, var = 'SampleID')

column_to_delete <- 1
pheno_2 <- pheno_2[, -column_to_delete, drop = FALSE]
head(pheno_2)

# check that cell types sum to 1
## use the dplyr:: because R was confused about pathways
cell_sum <- pheno_2 %>% 
  dplyr::select(c("cd8t", "cd4t", "nk", "bcell","mono", "neu")) %>% 
  rowSums()

## Get only the information we need
#Binary MDD state, sex, cell type proportion
cols_of_interest_ <- c("explainatory", "cd8t", "cd4t", "nk", "bcell","mono", "neu", "female")

pheno_2_C <- pheno_2[cols_of_interest]

prep pheno for current vs remitted}

## prep pheno for current vs remitted
## load phenotype data
pheno_3 <- read.csv('GTP_Current_Remitted_ASM_Pheno.csv', 
                    header = TRUE, 
                    sep = ",")

pheno_3 <- column_to_rownames(pheno_3, var = 'SampleID')

column_to_delete <- 1
pheno_3 <- pheno_3[, -column_to_delete, drop = FALSE]
head(pheno_3)

# check that cell types sum to 1
## use the dplyr:: because R was confused about pathways
cell_sum <- pheno_3 %>% 
  dplyr::select(c("cd8t", "cd4t", "nk", "bcell","mono", "neu")) %>% 
  rowSums()

## Get only the information we need
#Binary MDD state, sex, cell type proportion
cols_of_interest <- c("explainatory", "cd8t", "cd4t", "nk", "bcell","mono", "neu", "Sex")
pheno_3_C <- pheno_3[cols_of_interest]



#################################################### load DNHS CpGs for targeted replication}
## load DNHS files for targeted replication


cpg_sites_NC <- read.csv('DNHS_cpg_sites_list_pheno1.csv', 
                         header = TRUE, 
                         sep = ",")

cpg_sites_NR <- read.csv('DNHS_cpg_sites_list_pheno2.csv', 
                         header = TRUE, 
                         sep = ",")
cpg_sites_CR <- read.csv('DNHS_cpg_sites_list_pheno3.csv', 
                         header = TRUE, 
                         sep = ",")


DNHS_NC <- as.list(colnames(cpg_sites_NC))
GTP_NC <- rownames(myResults_sig_1)
NC_common <- intersect(DNHS_NC, GTP_NC)
print(NC_common)

DNHS_NR <- as.list(colnames(cpg_sites_NR))
GTP_NR <- rownames(myResults_sig_2)
NR_common <- intersect(DNHS_NR, GTP_NR)
print(NR_common)

DNHS_CR <- cpg_sites_CR$V1
GTP_CR <- rownames(myResults_sig_3)
CR_common <- intersect(DNHS_CR, GTP_CR)
print(CR_common)


#################################################### rank probes and run DMR analyses
## rank probes and run test never vs current

## Rank the probes in 2 groups as binary variable Yes/No
myRank_1 <- rankProbes(beta_1_top, 
                       pheno_1_C, 
                       #explanatory = 9, # this lets you pick the column that your explanatory var is in if it is not the first
                       refGroup = "control", 
                       typeInput = "beta")

## run the test and select the significant outputs

##safe <- gr_cpg_sites
myResults_1 <-run_mCSEA(myRank_1, beta_1_top,
                        pheno_1_C, cpg_site_list)

## found 50 DMRs with padj < 0.05

myResults_sig_1 <- get_significant(results = myResults_1)
dim(myResults_sig_1)
#dim(beta_wd_var)

View(myResults_sig_1)

## save results

write.csv(myResults_sig_1, file = "GTP_never_current_significant_mcSEA_Results", quote = FALSE)

## pull results/subset from myResults_sig from myResults so the genes can be plotted
gene_names_1 <- row.names(myResults_sig_1)
get_all_DMRs(myResults_1)

##adding chromosome info to results
ensembl <- useMart("ensembl")

library(org.Hs.eg.db)
library(Hmisc)
library(dbplyr)
# ensembl_host <- "https://www.ensembl.org/"

## Retrieve chromosome information for the gene names in myResults_sig

annots_1 <- select(org.Hs.eg.db, keys=rownames(myResults_sig_1), 
                   columns = c("SYMBOL", "GENENAME"), 
                   keytype = "SYMBOL")
print(annots_1)



## rank probes and run test never vs remitted

## Rank the probes in 2 groups as binary variable Yes/No
myRank_2 <- rankProbes(beta_2_top, 
                       pheno_2_C, 
                       #explanatory = 1, # this lets you pick the column that your explainatory var is in if it is not the first
                       refGroup = "control", 
                       typeInput = "beta")

## run the test and select the significant outputs
## First for all looking at DMR with treatment Treatment Response at timepoint 2 as a Binary variable
myResults_2 <-run_mCSEA(myRank_2, beta_2_top,
                        pheno_2_C, cpg_site_list)

## found 50 DMRs with padj < 0.05

myResults_sig_2 <- get_significant(results = myResults_2)
dim(myResults_sig_2)
#dim(beta_wd_var)

View(myResults_sig_2)

## save results

write.csv(myResults_sig_2, file = "GTP_never_remitted_significant_mcSEA_Results", quote = FALSE)

## pull results/subset from myResults_sig from myResults so the genes can be plotted
gene_names_2 <- row.names(myResults_sig_2)
get_all_DMRs(myResults_2)

## Retrieve chromosome information for the gene names in myResults_sig

annots_2 <- select(org.Hs.eg.db, keys=rownames(myResults_sig_2), 
                   columns = c("SYMBOL", "GENENAME"), 
                   keytype = "SYMBOL")
print(annots_2)



## rank probes and run test current vs remitted

## Rank the probes in 2 groups as binary variable Yes/No
myRank_3 <- rankProbes(beta_3_top, 
                       pheno_3_C, 
                       #explainatory = 62, # this lets you pick the column that your explainatory var is in if it is not the first
                       refGroup = "control", 
                       typeInput = "beta")

## run the test and select the significant outputs
## First for all looking at DMR with treatment Treatment Response at timepoint 2 as a Binary variable
myResults_3 <-run_mCSEA(myRank_3, beta_3_top,
                        pheno_3, cpg_site_list)


myResults_sig_3 <- get_significant(results = myResults_3)
dim(myResults_sig_3)
#dim(beta_wd_var)

View(myResults_sig_3)

## save results

write.csv(myResults_sig_3, file = "GTP_current_remitted_significant_mcSEA_Results", quote = FALSE)

## pull results/subset from myResults_sig from myResults so the genes can be plotted
gene_names_3 <- row.names(myResults_sig_3)
get_all_DMRs(myResults_3)

## Retrieve chromosome information for the gene names in myResults_sig

annots_3 <- select(org.Hs.eg.db, keys=rownames(myResults_sig_3), 
                   columns = c("SYMBOL", "GENENAME"), 
                   keytype = "SYMBOL")
print(annots_3)

################################################ print to table for word

## Install and load necessary packages
install.packages("tidyverse")
install.packages("flextable")
install.packages("officer")

library(tidyverse)
library(flextable)
library(officer)
library(dbplyr)

# Rename columns if necessary
#colnames(mcsea_results) <- c("GeneSet", "p_value", "score")

myResults_sig_1$RowNames <- rownames(myResults_sig_1)

# Reorder columns to have RowNames as the first column
myResults_sig_1 <- myResults_sig_1[, c(ncol(myResults_sig_1), 1:(ncol(myResults_sig_1)-1))]


# Create a flextable object
table <- flextable(myResults_sig_1)

# Customize table appearance if needed
table <- set_table_properties(table, layout = "autofit")

# Export the table to a Word document
doc <- read_docx()
doc <- doc %>% 
  body_add_flextable(value = table) %>%
  print(target = "mCSEA_Results_1.docx")




######################## Never Vs Remitted

myResults_sig_2$RowNames <- rownames(myResults_sig_2)

# Reorder columns to have RowNames as the first column
myResults_sig_2 <- myResults_sig_2[, c(ncol(myResults_sig_2), 1:(ncol(myResults_sig_2)-1))]


# Create a flextable object
table <- flextable(myResults_sig_2)

# Customize table appearance if needed
table <- set_table_properties(table, layout = "autofit")

# Export the table to a Word document
doc <- read_docx()
doc <- doc %>% 
  body_add_flextable(value = table) %>%
  print(target = "mCSEA_Results_2.docx")


######################## Current Vs Remitted

myResults_sig_3$RowNames <- rownames(myResults_sig_3)

# Reorder columns to have RowNames as the first column
myResults_sig_3 <- myResults_sig_3[, c(ncol(myResults_sig_3), 1:(ncol(myResults_sig_3)-1))]


# Create a flextable object
table <- flextable(myResults_sig_3)

# Customize table appearance if needed
table <- set_table_properties(table, layout = "autofit")

# Export the table to a Word document
doc <- read_docx()
doc <- doc %>% 
  body_add_flextable(value = table) %>%
  print(target = "mCSEA_Results_3.docx")

