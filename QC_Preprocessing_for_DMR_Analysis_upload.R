####### Quality control and preprocessing for DMR analyses

####################################### set working dirtectory and load packages

## set working directory
setwd('/Users/laniemullins/Desktop/MSPH_Thesis/Grady_Trauma_Project/')

## load packages
install.packages("illuminaio")
#remotes::install_github("illuminaio")
install.packages("https://github.com/hhhh5/ewastools/archive/master.tar.gz", repos = NULL, type = "source")
#devtools::install_github("hhhh5/ewastools")
library(ewastools)
install.packages("dplyr")
library(dplyr)
install.packages("data.table")
library(data.table)



######### Load Idats

# function to get idat files list
# Input: path to the directory having idats (subfolders of each chip)
# it will read names of all files, and create the paths to each file
# Output: list of paths

file_list <- function(path){
  files <- list.files(path, recursive = TRUE, pattern = ".idat")
  files <- gsub(paste(c("_Grn.idat", "_Red.idat"), collapse = "|"), "", files)
  file_paths <- paste0(path, files)
  files <- gsub(".*/","",files)
  return(list(names = unique(files), paths = unique(file_paths)))
}


############### Load sample sheet or phenotype file
# This file should include SentrixID and SentrixPosition as column names
# Then create SampleID by combining Sentrix barcode (Sentrix_ID) and Sentrix position (Sentrix_Position)
## did not need to do this for all of the samples, just the GTP 

pheno <- read.csv("GTP_ASM_Pheno.csv", header = TRUE, sep = ',', row.names = 1) # PHENOTYPE FILE - HERE

## already have SampleID do not need this step, but will leave just in case
#pheno$SampleID <- paste(pheno$SentrixID, #pheno$SentrixPosition,sep = "_") # create SampleID by combining Sentrix barcode and Sentrix position



######### read idats
## Path where you idat files are
## It can have many subfolders containing idats

## these files need to be unzipped, so in shell used:
##  gzip -d *.gz while in the desired directory

main_dir <- "/Users/laniemullins/Desktop/MSPH_Thesis/Grady_Trauma_Project/ASM_idats/"  # PATH TO IDATS FOLDER HERE, this is an example path, change accordingly

file_paths <- file_list(path = main_dir)

paths <- file_paths$paths[which(file_paths$names %in% pheno$SampleID)]

stopifnot(all(grepl(paste(c(pheno$SampleID), collapse = "|"), paths))) # stop if all samples are not in your path

meth <- read_idats(paths, quiet = FALSE)



######################################################### Control Metrics and Genotype/Sex Check

############################################# Control Metrics

## evaluates 17 control metrics which are describe in the BeadArray Controls Reporter Software Guide from Illumina.
## compares all 17 metrics against the thresholds recommended by Illumina.

ctrls <- control_metrics(meth)

pheno$failed = sample_failure(ctrls)
table(pheno$failed) # Samples with failed == TRUE should be identified and removed for further analysis
fails <- subset(pheno, failed == "TRUE")
message("Failed Samples: ", nrow(fails)) # No samples are failed


################# Genotype Calling and Outliers
# Genotype calling and outliers ----------------------------------------
# Need to do a quick pre-processing just for that step, but we won't be using these beta values for the further analysi


beta <- meth %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize
snps <- meth$manifest[probe_type == "rs", index]
snps <- beta[snps, ]

#' Check if you got only snps, if everything is fine, it should not give any error
stopifnot(all(grepl("rs", rownames(snps))))

# estimates the parameters of a mixture model consisting of three Beta distributions representing the heterozygous and the two homozygous genotypes,
# and a fourth component, a uniform distribution, representing outliers.
genotypes <- call_genotypes(snps, learn = FALSE)
# average log odds of belonging to the outlier component across all SNP probes.
# flagging samples with a score greater than -4 for exclusion is recommended, because of possible contamination.
pheno$outlier <- snp_outliers(genotypes)

# This will raise an error if 'pheno' is not a data table
# So we need to convert it to data table
if(!is.data.table(pheno)){
  pheno <- data.table(pheno)
}
pheno[outlier > -4, .(SampleID,outlier)]

# Flag the outlier samples: Flagged outlier samples are denoted as "Y", while good samples as "N"
pheno$outlierYN <- pheno$outlier > -4

# Check for duplicated samples
pheno$donor_id <- enumerate_sample_donors(genotypes)

# List duplicates
pheno[, n:=.N, by = donor_id]
pheno[n > 1, .(SampleID,donor_id)] # Check to see if there are any duplicates, and remove idats of the one from the folder

# Save the phenotype file with the QC details
write.csv(pheno, file = "GTP_ASM_Pheno_QC.csv", row.names = F)


load packages}

## clear room in environment 
rm(list=ls())
gc()

####################################### Start of Normalization and Preprocessing
###### install the required packages
bioc_packages <- c("BiocManager", "minfi", "IlluminaHumanMethylationEPICmanifest",
                   "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                   "sva", "limma", "impute", "BiocParallel")

cran_packages <- c( "data.table", "pamr", "feather", "tibble", "CpGassoc")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("Biomanager")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("sva")
BiocManager::install("limma")
BiocManager::install("impute")
install.packages("CpGassoc")
install.packages("data.table")
install.packages("pamr")
install.packages("tibble")
#install.packages("feather")


library(minfi)
library(CpGassoc)
library(impute)
library(pamr)
library(tibble)
library(data.table)
library(limma)

######################################################### Reading idats and loading needed QCd file
## function to read the idat files
# Input: path to the idats, Phenotype file
# Output: RGset, Mset, GMset in a list
read_idat_files <- function(path, pheno_file){
  RGset <- read.metharray.exp(path, recursive = TRUE, verbose = TRUE, force = TRUE)
  RGset <- RGset[, colnames(RGset) %in% pheno$SampleID]
  
  if(!all(colnames(RGset) %in% pheno$SampleID)){
    stop("samples in idats and sample sheet are not same")
  }
  
  Mset <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod="single")
  GMset <- mapToGenome(Mset, mergeManifest = FALSE)
  
  return(list(RGset = RGset, Mset =  Mset, GMset = GMset))
}

## Load the new phenotype file with QC info
## Remove failed samples
pheno <- read.csv("GTP_ASM_Pheno_QC.csv",as.is = T) # WRITE PHENOTYPE FILE HERE
pheno <- subset(pheno, failed == "FALSE")

## This is the path to the folders having idats
main_dir <- "/Users/laniemullins/Desktop/MSPH_Thesis/Grady_Trauma_Project/ASM_idats/"  # PATH TO IDATS FOLDER HERE, this is an example path, change accordingly  

## Now call the function
## you will get RGset, Mset and GMset in a list
output <- read_idat_files(path = main_dir, pheno_file = pheno)

############################################## Sex Check and Sample Preprocessing

## predict the sex
## function to match sex with predicted sex
## remove the sex mismatching samples from phenotype file, RGset, Mset and GMset
## Input: list of (RGset, Mset and GMset), phenotype file, name of sex column in your phenotype file
## Sex should be coded as follows: 0 or M for males and 1 or F for females  
## Output : predicted sex and sex mismatches

check_sex_info <- function(inputdata , pheno_file, pheno_sex_col){
  message("Processing, please wait ...")
  pheno_file[,pheno_sex_col][pheno_file[,pheno_sex_col] == 0] <- "M"
  pheno_file[,pheno_sex_col][pheno_file[,pheno_sex_col] == 1] <- "F"
  predicted_sex <- getSex(object = inputdata$GMset, cutoff = -2)
  
  rownames(pheno_file) <- pheno_file$SampleID
  pheno <- pheno_file[order(rownames(pheno_file)), ]
  sex <- predicted_sex[order(rownames(predicted_sex)), ]
  
  if(!all(rownames(pheno) == rownames(sex))){
    stop("Phenotype and predicted sex samples doesn't match")
  }
  RGset <- inputdata$RGset
  Mset <- inputdata$Mset
  GMset <- inputdata$GMset
  mm <- rownames(pheno[pheno[[pheno_sex_col]] != sex$predictedSex, ])
  message("Samples with sex mismatch: ", length(mm))
  if(length(mm)){
    RGset <- RGset[, which(!colnames(RGset) %in% mm)]
    Mset <- Mset[, which(!colnames(Mset) %in% mm)]
    GMset <- GMset[, which(!colnames(GMset) %in% mm)]
    pheno <- pheno[which(!pheno$SampleID %in% mm), ]
  }
  return(list(pheno = pheno, RGset = RGset, Mset = Mset, GMset = GMset, predicted_sex = predicted_sex))
}

## call the function to check sex information
## pheno_sex_col = the name of the sex column in the phenotype (Sex, Gender etc)
all_info <- check_sex_info(inputdata = output, pheno_file = pheno, pheno_sex_col = "Sex")

## check if all data have same samples
## all should be true
lapply(all_info[c(2:4)], function(x) all(colnames(x) %in% all_info$pheno$SampleID))

## write predicted sex onto the phenotype file and flag the mismatches
pheno <- pheno[order(pheno$SampleID),]
all_info$predicted_sex <- all_info$predicted_sex[order(rownames(all_info$predicted_sex)),]
all(rownames(all_info$predicted_sex)==pheno$SampleID) # This should be TRUE
pheno$predictedSex <- all_info$predicted_sex$predictedSex
if(all(0|1%in%pheno$Sex) == T) {
  pheno$predictedSex[pheno$predictedSex == "M"] <- 0 
  pheno$predictedSex[pheno$predictedSex == "F"] <- 1
}
pheno$SexMismatch <- pheno$Sex != pheno$predictedSex 

## extract preprocessed pval, beta, methylated and unmethylated signal files

pval <- detectionP(all_info$RGset)
beta <- getBeta(all_info$Mset)
signalA <- getUnmeth(all_info$Mset)
signalB <- getMeth(all_info$Mset)

## save preprocessed pval, beta, methylated and unmethylated signal files
save(pval, beta, signalA, signalB, file = "IntermediaryFiles.RData")

## save the new phenotype file with flagged sex mismatches
write.csv(pheno, file = "GTP_ASM_Pheno_QC_postMINFI.csv", row.names = F)

## remove all objects from your workspace
rm(list=ls())
gc()

######################################################### QC with CpGassoc

## load the preprocessed pval, beta, methylated and unmethylated signal files

load("IntermediaryFiles.RData")

beta.new <- cpg.qc(beta, signalA, signalB, pval, p.cutoff=.01,cpg.miss=.1, sample.miss=.1)

#[1] "Removed x samples with low signal"
#[1] "Removed x CpG sites with missing data for > 0.1 of samples"
#[1] "Removed x samples with missing data for > 0.1 of CpG sites"

## Remove cross reactive probes: Got the cross reactive probes from a paper (http://www.sciencedirect.com/science/article/pii/S221359601630071X).

## processing QCd data to remove only Cross Reactive Probes
cross <- read.csv("cross.csv", stringsAsFactors = FALSE, header = TRUE) # This file is available with the pipeline package we've shared
cross <- cross[, 1]
beta.new <- beta.new[!(row.names(beta.new) %in% cross), ]

## run test if all cross cpgs are removed
stopifnot(all(!cross %in% rownames(beta.new)))

## save Noob Normalized Cross-Reactive Probes removed beta values
write.csv(beta.new, file = "GTP_ASM_noob_qcd_crossReactiveProbesRemoved.csv", quote = FALSE, row.names = TRUE) 
# keep the noob normalized and QCd data

#################################### ComBAT normalization to adjust for batch effects (chip and array position

## clean/clear up space
rm(list=ls())
gc()

## reinstall packages if necessary
install.packages("data.table")
library(data.table)

## load Noob Normalized Cross-Reactive Probes removed beta values
beta <- fread("GTP_ASM_noob_qcd_crossReactiveProbesRemoved.csv", data.table = F) # cpgs x samples
row.names(beta)<-beta$V1
beta<-beta[,-1]

## load the new phenotype file with QC info
phen <- read.csv("GTP_ASM_Pheno_QC_postMINFI.csv", stringsAsFactors = FALSE, header = TRUE)
row.names(phen)<-phen$SampleID

## remove the failed samples and sex mismatches
phen <- subset(phen, failed == "FALSE" & SexMismatch == "FALSE")
shuffle <- sample(nrow(phen))
phen <- phen[shuffle, ]

## mthylation IDs (Row names of phenotype file and Column names of beta matrix) should match
beta <- beta[, which(names(beta) %in% row.names(phen))]
phen <- phen[row.names(phen) %in% names(beta), ]
beta <- beta[, order(names(beta))]
phen <- phen[order(row.names(phen)), ]
stopifnot(all(colnames(beta) == rownames(phen)))

## Define Variables for Combat step
## Define model.matrix, which includes the variables that you want to protect when adjusting for chip and position
## Generally variables that we use as covariates in the EWAS (sex, age, main phenotype -PTSD-, smoking) are included in the model.matrix
#sex <- "Sex"         # Name of Sex Variable: Males coded as 0, females coded as 1
## sex removed from GTP analysis because it confounded batch effects
age <- "Age"         # Name of Age Variable
race <- "race_ethnic"
PTSD <- "PTSD"       # Name of PTSD Variable: Cases coded as 1, controls coded as 0
chip <- "sentrixID"
position <- "sentrixPosition"

## You should not have NAs in model matrix, so we remove subjects with no phenotype info
print(paste0("Samples with no PTSD information = ", sum(is.na(phen[,PTSD]))))
#print(paste0("Samples with no Sex information = ", sum(is.na(phen[,sex]))))
print(paste0("Samples with no Age information = ", sum(is.na(phen[,age]))))
print(paste0("Samples with no Race information = ", sum(is.na(phen[,race]))))


naindex <- (!is.na(phen[, race]) & !is.na(phen[, age]) & !is.na(phen[, PTSD]))
phen <- phen[naindex, ]
beta <- beta[, naindex]
stopifnot(all(colnames(beta) == rownames(phen)))

chip <- as.factor(phen[,chip])
position <- as.factor(phen[,position])
ptsd<-as.factor(phen[,PTSD])
age<-as.numeric(phen[,age])
#sex<-as.factor(phen[,sex])
race<-as.factor(phen[,race])


moddata <- model.matrix(~age+race+PTSD) 

## remaining samples
print(paste0("Remaining Samples = ", nrow(phen))) # N = 142,  2 failed QC

## ComBAT does not handle NAs in the methylation file,
## so we have to impute NAs in the methylation beta matrix
## Log transform to normalize the data (we'll reverse that later)
beta <- log((beta/(1-beta)))
beta.imputed <- impute.knn(as.matrix(beta))
beta.imputed <- beta.imputed$data

## clear space before running the next step if needed
rm(beta)
gc()


## load packages needed because of package nesting
install.packages("BiocParallel")
library(BiocParallel)

# Run ComBat
install.packages("sva")
library(sva)
SerialParam()
combat_beta <- ComBat(dat = beta.imputed, mod = moddata, batch = chip, BPPARAM = SerialParam())
combat_beta <- ComBat(dat = combat_beta, mod = moddata, batch = position, BPPARAM = SerialParam())
## found 85 batches


## reverse Beta Values
reversbeta <- 1/(1+(1/exp(combat_beta)))

## need to put NAs back (from the original matrix) to ComBAT adjusted beta matrix
## don't want to use imputed beta values for missing data
norm <- fread("GTP_ASM_noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(norm) <- norm$V1
norm<-norm[,-1]

beta <- as.data.frame(reversbeta)

norm <- norm[, names(norm) %in% names(beta)]
norm <- norm[, order(names(norm))]
beta <- beta[, order(names(beta))]
table(names(norm) == names(beta))

beta[is.na(norm)] <- NA
beta[beta >= 0.9999998] <- 0.9999999

write.csv(beta,file="GTP_ASM_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_race.csv",quote=FALSE,row.names=TRUE)


## cell proportion was already calculated and within pheno file on GEO so it was not included here
