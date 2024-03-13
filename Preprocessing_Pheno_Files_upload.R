######## preprocessing pheno files 

#### Pulling and Splitting Participants by MDD Status + Age and Sex Matching

########################################################### set WD and load files

setwd('/Users/laniemullins/Desktop/School/MSPH_Genomics/Uddin Lab /DNHS_Depression/')

## read csv with covariates

methylation_data <- read.csv('DNHS_QCd_TraumaExp_CovarAdjusted.csv', header = TRUE, sep = ",")


## matching participants with EPIC data to MDD binary data
# load data
unique_participants_EPIC <- read.csv('DNHS_Epic_Phenotypes_Unique_Ppts.csv', 
                                     header = TRUE, 
                                     sep = ",")
dnhs_MDD_data <- read.csv('DNHS_MDD.csv', 
                          header = TRUE, 
                          sep = ',')


########################################################### find unique participants

## find intersect
RESP_MDD <- dnhs_MDD_data$RESP
RESP_EPIC <- unique_participants_EPIC$RESP
commonRESP <- intersect(RESP_EPIC, RESP_MDD)

## in a new data frame, filter the MDD data by the common RESP IDS.
RESP_in_EPICandMDD <- dnhs_MDD_data[dnhs_MDD_data$RESP %in% commonRESP, ]

## looking for RESP, current, and lifetime MDD (columns, 1-3, 10, and 17)
## this is all of the data including repeated RESP for the EPIC/MDD

dnhs_MDD_ROI <- RESP_in_EPICandMDD[, c(1:3, 10, 17)]

dnhs_MDD_ROI_file <- write.csv(dnhs_MDD_ROI, file = 'dnhsMDD_data_ROI.csv', quote = F)
dnhs_MDD_ROI_file <- read.csv('dnhsMDD_data_ROI.csv', 
                              header = TRUE,
                              row.name = 1,
                              sep = ',')
uniqueParticipants <- dnhs_MDD_ROI_file[!duplicated(dnhs_MDD_ROI_file$RESP), ]
## nrow(uniqueParticipants = 530)

## now we have a file that has only the RESP IDs that the EPIC and MDD file have in common

seperate into MDD disease states}

####################################### REMITTED START
## filter out individuals that are remitted (lifetime is a yes, and current is a no)
dnhs_MDD_Remitted <- subset(dnhs_MDD_ROI_file, (dnhs_MDD_ROI_file$LifetimeMDD == 1 & dnhs_MDD_ROI_file$CurrentMDD == 0))
head(dnhs_MDD_Remitted)
dnhs_MDD_Remitted_File <- write.csv(dnhs_MDD_Remitted, file = 'dnhs_MDD_remitted.csv', quote = F)
dnhs_MDD_Remitted_File <- read.csv('dnhs_MDD_remitted.csv', 
                                   header = TRUE,
                                   row.name = 1,
                                   sep = ',')

# find the unique RESP for Remitted
dnhs_MDD_Remitted_Unique <- dnhs_MDD_Remitted_File[!duplicated(dnhs_MDD_Remitted_File$RESP), ]
dnhs_MDD_Remitted_Unique <- write.csv(dnhs_MDD_Remitted_Unique, file = 'dnhs_Remitted_Unique.csv', quote = F)
dnhs_MDD_Remitted_Unique <- read.csv('dnhs_Remitted_Unique.csv', 
                                     header = TRUE,
                                     row.name = 1,
                                     sep = ',')
head(dnhs_MDD_Remitted_Unique)
nrow(dnhs_MDD_Remitted_Unique)
## Remitted (n = 80)


############################################# CURRENT START
dnhs_MDD_Current <- subset(dnhs_MDD_ROI_file, (dnhs_MDD_ROI_file$LifetimeMDD ==1 & dnhs_MDD_ROI_file$CurrentMDD == 1))
dnhs_MDD_Current <- write.csv(dnhs_MDD_Current, file = 'dnhs_MDD_current.csv', quote = F)
dnhs_MDD_Current <- read.csv('dnhs_MDD_current.csv',
                             header = TRUE,
                             row.name = 1,
                             sep = ',')
head(dnhs_MDD_Current)
dnhs_MDD_Current_Unique <- dnhs_MDD_Current[!duplicated(dnhs_MDD_Current$RESP), ]
dnhs_MDD_Current_Unique <- write.csv(dnhs_MDD_Current_Unique, file = 'dnhs_Current_Unique.csv', quote = F)
dnhs_MDD_Current_Unique <- read.csv('dnhs_Current_Unique.csv', 
                                    header = TRUE,
                                    row.name = 1,
                                    sep = ',')
head(dnhs_MDD_Current_Unique)
nrow(dnhs_MDD_Current_Unique)
## Current (n = 49)

############################################## LIFETIME NEVER START
dnhs_MDD_Never <- subset(dnhs_MDD_ROI_file, (dnhs_MDD_ROI_file$LifetimeMDD == 0 & dnhs_MDD_ROI_file$CurrentMDD == 0))
dnhs_MDD_Never <- write.csv(dnhs_MDD_Never, file = 'dnhs_MDD_Never.csv', quote = F)
dnhs_MDD_Never <- read.csv('dnhs_MDD_Never.csv',
                           header = TRUE,
                           row.name = 1,
                           sep = ',')
head(dnhs_MDD_Never)
dnhs_MDD_Never_Unique <- dnhs_MDD_Never[!duplicated(dnhs_MDD_Never$RESP), ]
dnhs_MDD_Never_Unique <- write.csv(dnhs_MDD_Never_Unique, file = 'dnhs_MDD_Never_Unique.csv', quote = F)
dnhs_MDD_Never_Unique <- read.csv('dnhs_MDD_Never_Unique.csv', 
                                  header = TRUE,
                                  row.name = 1,
                                  sep = ',')
head(dnhs_MDD_Never_Unique)
nrow(dnhs_MDD_Never_Unique)
## Never (n = 461)

#### Summary: Total (n=598), Never (n = 461), Remitted = (n = 80), Current (n = 49)
## Note to self, there may be individuals over lapping in these cohorts, check and favor keeping current/remitted

## download package to use antijoin
library(dplyr)

#################################### CHECKING FOR INTERSECT
common_entries_R_C <- intersect(dnhs_MDD_Current_Unique$RESP, dnhs_MDD_Remitted_Unique$RESP)
common_entries_R_C <- write.csv(common_entries_R_C, file = 'common_entries_remitted_current.csv',quote = F)
common_entries_R_C <- read.csv('common_entries_remitted_current.csv', header = TRUE, row.names = 1, sep = ',')
colnames(common_entries_R_C) <- c("RESP")
nrow(common_entries_R_C)

## NOTE: remove from which for current and remitted (find intersection with remitted and remove them from the file, they are stored in the csv created above)

ids_to_remove_R_C <- data.frame(common_entries_R_C$RESP)
colnames(ids_to_remove_R_C) <- c('RESP')
filtered_remitted_data <- dnhs_MDD_Remitted_Unique[!dnhs_MDD_Remitted_Unique$RESP %in% common_entries_R_C$RESP, ]
nrow(filtered_remitted_data)
filtered_remitted_data <- write.csv(filtered_remitted_data, file = 'filtered_remitted_data.csv', quote = F)

## - 25 from remitted (n = 55)
## unique participants with filtered out commonalities are in "filtered_remitted_data"


########### next set of overlap checking
common_entries_C_N <- intersect(dnhs_MDD_Current_Unique$RESP, dnhs_MDD_Never_Unique$RESP)
common_entries_C_N <- write.csv(common_entries_C_N, file = 'common_entries_current_never.csv',quote = F)
common_entries_C_N <- read.csv('common_entries_current_never.csv', header = TRUE, row.names = 1, sep = ',')
colnames(common_entries_C_N) <- c("RESP")
nrow(common_entries_C_N)


ids_to_remove_C_N <- data.frame(common_entries_C_N$RESP)
colnames(ids_to_remove_C_N) <- c('RESP')
filtered_never_data <- dnhs_MDD_Never_Unique[!dnhs_MDD_Never_Unique$RESP %in% common_entries_C_N$RESP, ]
nrow(filtered_never_data)
#filtered_never_data <- write.csv(filtered_never_data, file = 'filtered_never_data.csv', quote = F)

## see note above
## - 21 from never (n = 440)

common_entries_R_N <- intersect(dnhs_MDD_Remitted_Unique$RESP, dnhs_MDD_Never_Unique$RESP)
common_entries_R_N <- write.csv(common_entries_R_N, file = 'common_entries_remitted_never.csv',quote = F)
common_entries_R_N <- read.csv('common_entries_remitted_never.csv', header = TRUE, row.names = 1, sep = ',')
colnames(common_entries_R_N) <- c('RESP')
nrow(common_entries_R_N)


ids_to_remove_R_N <- data.frame(common_entries_R_N$RESP)
colnames(ids_to_remove_R_N) <- c('RESP')
filtered_never_data2 <- dnhs_MDD_Never_Unique[!dnhs_MDD_Never_Unique$RESP %in% common_entries_R_N$RESP, ]
nrow(filtered_never_data2)

## compare the 2 filtered Never data sets and only keep the ones that DO intersect 
## the ones that do NOT intersect overlapped with either C or R and were removed during filtering
ids_to_remove_N <- data.frame(filtered_never_data2$RESP)
colnames(ids_to_remove_N) <- c("RESP")
nrow(ids_to_remove_N)


filtered_never_data <- filtered_never_data[filtered_never_data$RESP %in% ids_to_remove_N$RESP, ]
nrow(filtered_never_data)
filtered_never_data <- write.csv(filtered_never_data, file = 'filtered_never_data.csv', quote = F)
## (n = 426)

## find the intersect of the common entry files and the unique participant files and remove as described
## check to make sure there is no longer any intersection between the files for Unique participants

#intersect_check1 <- intersect(dnhs_MDD_Current_Unique$RESP, filtered_never_data$RESP)
#nrow(intersect_check1)


## Note: all 3 possible intersections were checked, and there was no overlap in any of the groups

## files of interest are read below
filtered_never_data <- read.csv('filtered_never_data.csv', header = TRUE, row.names = 1, sep= ',')
filtered_remitted_data <- read.csv('filtered_remitted_data.csv', header = TRUE, row.names = 1, sep = ',')
filtered_current_data <- read.csv('dnhs_Current_Unique.csv', header = TRUE, row.names = 1, sep = ',')

## Summary: Total Unique Participants (n = 530), Never (n = 426), Remitted (n = 55), Current (n = 49)


################################################# age+sex matching

library(dplyr)

## age and sex matching to keep sample sizes consistent
merged_never_remitted <- merge(filtered_never_dataDone, filtered_remitted_dataDone, by = c("female", "Age"))

## delete rows that have a duplicated RESP.y (remitted represented more than once)

duplicated_rows_RN_R <- duplicated(merged_never_remitted$RESP.y)
merged_never_remitted_no_dups <- merged_never_remitted[!duplicated_rows_RN_R, ]

## delete rows that have a duplicated RESP.x(never represented more than once)
duplicated_rows_RN_N <- duplicated(merged_never_remitted_no_dups$RESP.x)
merged_never_remitted_no_dups_2 <- merged_never_remitted_no_dups[!duplicated_rows_RN_N, ]
## n = 41, so there are exclusions. 

## get rid of all of the columns for remitted to leave me with a never subset
columns_to_remove_never_A <- c("X.y", "BloodID.y", "RESP.y", "wave.y", "LifetimeMDD.y", "CurrentMDD.y", "SampleID.y", "CD8T.y", "CD4T.y","NK.y", "Bcell.y", "Mono.y", "Neu.y", "Comp.2.y", "Comp.3.y", "SmoS.y")
never_age_sex_matched_remitted <- merged_never_remitted_no_dups_2[, !colnames(merged_never_remitted_no_dups_2) %in% columns_to_remove_never_A]

columns_to_remove_remitted <- c("X.x", "BloodID.x", "RESP.x", "wave.x", "LifetimeMDD.x", "CurrentMDD.x", "SampleID.x", "CD8T.x", "CD4T.x","NK.x", "Bcell.x", "Mono.x", "Neu.x", "Comp.2.x", "Comp.3.x", "SmoS.x")
remitted_age_sex_matched <- merged_never_remitted_no_dups_2[, !colnames(merged_never_remitted_no_dups_2) %in% columns_to_remove_remitted]

file_never_ageSex_Matched_remitted <- write.csv(never_age_sex_matched_remitted, file = 'Never_Remitted_Age_Sex_Match.csv', quote = F)

file_remitted_ageSex_Matched <- write.csv(remitted_age_sex_matched, file = 'Remitted_Age_Sex_Match.csv', quote = F)


## remove the never_remitted matches from 
never_filtered_1 <- subset(filtered_never_dataDone, !(RESP %in% never_age_sex_matched_remitted$RESP.x))
merge_never_current <- merge(never_filtered_1, filtered_current_dataDone, by =c("female", "Age"))

## delete rows that have a duplicated RESP.y (remitted represented more than once)

duplicated_rows_CN_C <- duplicated(merge_never_current$RESP.y)
merged_never_current_no_dups <- merge_never_current[!duplicated_rows_CN_C, ]

## delete rows that have a duplicated RESP.x(never represented more than once)
duplicated_rows_CN_N <- duplicated(merged_never_current_no_dups$RESP.x)
merged_never_current_no_dups_2 <- merged_never_current_no_dups[!duplicated_rows_CN_N, ]

## n = 31, so there are exclusions. 
## get rid of all of the columns for current to leave me with a never subset

columns_to_remove_never_B <- c("X.y", "BloodID.y", "RESP.y", "wave.y", "LifetimeMDD.y", "CurrentMDD.y", "SampleID.y", "CD8T.y", "CD4T.y","NK.y", "Bcell.y", "Mono.y", "Neu.y", "Comp.2.y", "Comp.3.y", "SmoS.y")

never_age_sex_matched_current <- merged_never_current_no_dups_2[, !colnames(merged_never_current_no_dups_2) %in% columns_to_remove_never_B]

columns_to_remove_current <- c("X.x", "BloodID.x", "RESP.x", "wave.x", "LifetimeMDD.x", "CurrentMDD.x", "SampleID.x", "CD8T.x", "CD4T.x","NK.x", "Bcell.x", "Mono.x", "Neu.x", "Comp.2.x", "Comp.3.x", "SmoS.x")

current_age_sex_matched <- merged_never_current_no_dups_2[, !colnames(merged_never_current_no_dups_2) %in% columns_to_remove_current]

file_never_ageSex_Matched_current <- write.csv(never_age_sex_matched_current, file = 'Never_Current_Age_Sex_Match.csv', quote = F)

file_current_ageSex_Matched <- write.csv(current_age_sex_matched, file = 'Current_Age_Sex_Match.csv', quote = F)

all_never <- rbind(never_age_sex_matched_current, never_age_sex_matched_remitted)

file_all_all_never <- write.csv(all_never, file = 'Final_Never_Age_Sex_Match.csv', quote = F)

## now the total number of participants is 144
##### need these two files to be separate to find the corresponding beta values for the analysis

################################################################ split methylation file

## load needed package
library(dplyr)

## load file of beta values
beta_Values <- read.csv('DNHS_QCd_TraumaExp_CovarAdjusted.csv', 
                        header = TRUE,
                        row.names = 1,
                        sep = ",")

## there is an X in front of all of my column names for some reason, so lets remove them

colnames(beta_Values)
colnames(beta_Values) <- substring(colnames(beta_Values), 2)
colnames(beta_Values)
head(beta_Values)
rownames(beta_Values)


casm <- read.csv('Current_Age_Sex_Match.csv', 
                 header = TRUE, 
                 sep = ",")

rasm <- read.csv('Remitted_Age_Sex_Match.csv', 
                 header = TRUE, 
                 sep = ",")
nasm_1 <- read.csv('Never_Current_Age_Sex_Match.csv', 
                   header = TRUE, 
                   sep = ",")

nasm_2 <- read.csv('Never_Remitted_Age_Sex_Match.csv', 
                   header = TRUE, 
                   sep = ",")



############### get desired columns CURRENT

# for some reason, this removed the first column which is the name of my CpG sites...
# subset the columns of the second file to remove the unwanted columns
colnames(casm)[colnames(casm) == "SampleID.y"] <- "SampleID"
current_samples_casm <- casm$SampleID
current_MDD_beta_values_casm <- beta_Values[ , current_samples_casm] ## indexing fixed all of my problems
casm_current_MDD_beta_values_FILE <- write.csv(current_MDD_beta_values_casm, file = 'current_MDD_beta_values_Age_Sex_Match.csv', quote = F)

################# get desired columns REMITTED

colnames(rasm)[colnames(rasm) == "SampleID.y"] <- "SampleID"
remitted_samples_rasm <- rasm$SampleID
remitted_MDD_beta_values_rasm <- beta_Values[ , remitted_samples_rasm] 
rasm_remitted_MDD_beta_values_FILE <- write.csv(remitted_MDD_beta_values_rasm, file = 'remitted_MDD_beta_values_Age_Sex_Match.csv', quote = F)

################# get desired columns NEVER

## doing this twice because never is ASM twice

colnames(nasm_1)[colnames(nasm_1) == "SampleID.y"] <- "SampleID"
never_samples_nasm_1 <- nasm_1$SampleID
never_MDD_beta_values_nasm_1 <- beta_Values[ , never_samples_nasm_1] 
nasm_1_never_MDD_beta_values_FILE <- write.csv(never_MDD_beta_values_nasm_1, file = 'never_MDD_beta_values_Age_Sex_Match_to_current.csv', quote = F)


colnames(nasm_2)[colnames(nasm_2) == "SampleID.y"] <- "SampleID"
never_samples_nasm_2 <- nasm_2$SampleID
never_MDD_beta_values_nasm_2 <- beta_Values[ , never_samples_nasm_2] 
nasm_2_never_MDD_beta_values_FILE <- write.csv(never_MDD_beta_values_nasm_2, file = 'never_MDD_beta_values_Age_Sex_Match_to_remitted.csv', quote = F)

####### NOTE:
## may need to further clean/ join these

########################################################### joining and cleaning file formatting

## read files
casm <- read.csv('Current_Age_Sex_Match.csv', 
                 header = TRUE, 
                 sep = ",")

rasm <- read.csv('Remitted_Age_Sex_Match.csv', 
                 header = TRUE, 
                 sep = ",")
rasm_2 <- read.csv('Remitted_Age_Sex_Match_2.csv', 
                   header = TRUE, 
                   sep = ",")
nasm_1 <- read.csv('Never_Current_Age_Sex_Match.csv', 
                   header = TRUE, 
                   sep = ",")
nasm_2 <- read.csv('Never_Remitted_Age_Sex_Match.csv', 
                   header = TRUE, 
                   sep = ",")

beta_rasm <- read.csv('remitted_MDD_beta_values_Age_Sex_Match.csv', 
                      header = TRUE, 
                      sep = ",")

beta_casm <- read.csv('current_MDD_beta_values_Age_Sex_Match.csv', 
                      header = TRUE, 
                      sep = ",")

beta_nasm_1 <- read.csv('never_MDD_beta_values_Age_Sex_Match_to_current.csv', 
                        header = TRUE, 
                        sep = ",")

beta_nasm_2 <- read.csv('never_MDD_beta_values_Age_Sex_Match_to_remitted.csv', 
                        header = TRUE, 
                        sep = ",")


## define new columns
case <- c("case")
control <- c("control")

## add columns indicating if sample is case or control
library(dplyr)

casm <- casm %>% select(RESP.y, everything())
print(casm)
casm <- cbind(case, casm)
colnames(casm)[1] <- "explainatory"

nasm_1 <- nasm_1 %>% select(RESP.x, everything())
colnames(nasm_1)[1] <- "RESP.y"
print(nasm_1)
nasm_1 <- cbind(control, nasm_1)
print(nasm_1)
colnames(nasm_1)[1] <- "explainatory"

nasm_2 <- nasm_2 %>% select(RESP.x, everything())
colnames(nasm_2)[1] <- "RESP.y"
print(nasm_2)
nasm_2 <- cbind(control, nasm_2)
print(nasm_2)
colnames(nasm_2)[1] <- "explainatory"

rasm <- rasm %>% select(RESP.y, everything())
print(rasm)
rasm_case <- cbind(case, rasm)
print(rasm_case)
colnames(rasm_case)[1] <- "explainatory"

rasm_2 <- rasm_2 %>% select(RESP.y, everything())
print(rasm_2)
rasm_control <- cbind(control, rasm_2)
print(rasm_control)
colnames(rasm_control)[1] <- "explainatory"

## clean names for merging
# casm
column_to_delete <- "X.y"
casm <- casm[, !(colnames(casm) %in% column_to_delete), drop = FALSE]
column_to_delete <- "X"
casm <- casm[, !(colnames(casm) %in% column_to_delete), drop = FALSE]
new_column_names <- c("explainatory", "RESP", "female", 
                      "Age", "BloodID", "wave", "LifetimeMDD", "CurrentMDD",
                      "SampleID", "CD8T", "CD4T", "NK",
                      "Bcell", "Mono", "Neu", "Comp.2",
                      "Comp.3", "SmoS")
colnames(casm) <- new_column_names

# nasm_1
column_to_delete <- "X.x"
nasm_1 <- nasm_1[, !(colnames(nasm_1) %in% column_to_delete), drop = FALSE]
column_to_delete <- "X"
nasm_1 <- nasm_1[, !(colnames(nasm_1) %in% column_to_delete), drop = FALSE]
colnames(nasm_1) <- new_column_names
head(nasm_1)

# nasm_2
column_to_delete <- "X.x"
nasm_2 <- nasm_2[, !(colnames(nasm_2) %in% column_to_delete), drop = FALSE]
column_to_delete <- "X"
nasm_2 <- nasm_2[, !(colnames(nasm_2) %in% column_to_delete), drop = FALSE]
colnames(nasm_2) <- new_column_names
head(nasm_2)

# rasm_case
column_to_delete <- "X.y"
rasm_case <- rasm_case[, !(colnames(rasm_case) %in% column_to_delete), drop = FALSE]
column_to_delete <- "X"
rasm_case <- rasm_case[, !(colnames(rasm_case) %in% column_to_delete), drop = FALSE]
colnames(rasm_case) <- new_column_names

# rasm_control
column_to_delete <- "X.y"
rasm_control <- rasm_control[, !(colnames(rasm_control) %in% column_to_delete), drop = FALSE]
column_to_delete <- "X"
rasm_control <- rasm_control[, !(colnames(rasm_control) %in% column_to_delete), drop = FALSE]
colnames(rasm_control) <- new_column_names


## merge files of interest into 3 pheno and beta files 
## pheno_1/beta_1 is never and current, pheno_2/beta_2 is never, remitted, pheno_3/beta_3 is remitted, current

pheno_1 <- rbind(nasm_1, casm)
pheno_1_FILE <- write.csv(pheno_1, file = 'pheno_1.csv', quote = F)

pheno_2 <- rbind(nasm_2, rasm_case)
pheno_2_FILE <- write.csv(pheno_2, file = 'pheno_2.csv', quote = F)

pheno_3 <- rbind(rasm_control, casm)
pheno_3_FILE <- write.csv(pheno_3, file = 'pheno_3.csv', quote = F)

## merge beta_value files
## Note: cbind gives me the weird X, must fix before writing files

beta_1 <- cbind(beta_nasm_1, beta_casm)
colnames(beta_1) <- substring(colnames(beta_1), 2)
head(beta_1)
beta_1_FILE <- write.csv(beta_1, file = 'beta_1.csv', quote = F)


beta_2<- cbind(beta_nasm_2, beta_rasm)
colnames(beta_2) <- substring(colnames(beta_2), 2)
beta_2_FILE <- write.csv(beta_2, file = 'beta_2.csv', quote = F)

beta_3 <- cbind(beta_rasm, beta_casm)
colnames(beta_3) <- substring(colnames(beta_3), 2)
beta_3_FILE <- write.csv(beta_3, file = 'beta_3.csv', quote = F)


