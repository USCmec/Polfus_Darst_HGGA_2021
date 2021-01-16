##############################################################################################################################################
# This is example code for the PRS analysis found in Polfus, Darst et al., HGG Advances
# 
# This code provides an example of how the multiethnic PRS was calculated and how association tests were carried out within each population.
# Input files are described but cannot be provided due to IRB restrictions on individual level data.
#
# Author: Burcu F. Darst
# Date: 01/16/2021


### Load R libraries
library(data.table)


### Load genotype and phenotype files
## Participants must be *ordered identically* between the two files; this code assumes each file is limited to one population
# Genotype file is coded to reflect the number of risk-increasing alleles carried by each participant
# The variants included in this file and the risk-increasing allele of each variant are based on Vujkovic et al., Nature Genetics 2020
# This is an NxM genotype file, where N are the participants (rows) and M are the variants (columns)
gen = fread("genotypeFile.csv")
# This is an NxP genotype file, where N are the partcipants (rows) and P are the phenotypes (columns)
phen = fread("phenotypeFile.txt")


### Load summary statistics
# This is a vector of variant-specific weights, where variants are *ordered identically* as the "gen" file
# Weights for this analysis were obtained from Vujkovic et al., Nature Genetics 2020
# Weights are all coded to reflect the risk-increasing allele
weights = fread("weights.csv")


### Calculate PRS (weighted sum of risk variants, where weights are variant-specific)
gen = as.matrix(gen)
PRS = gen%*%weights


### Convert PRS to categorical variable
PRS.categories <-c(0, .1, .2, .3, .4, .6, .7, .8, .9,  1)
PRS.ref <- "40% - 60%"
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
# Decile cutoffs are based on controls within each population
q.v <- quantile(PRS[phen$Y==0], probs=PRS.categories)
q.v[1] <- min(q.v[1], min(PRS)) # if min(PRS) not in controls
q.v[length(q.v)] <- max(q.v[length(q.v)], max(PRS)) # if max(PRS) is not in controls
PRS.c <- cut(PRS, breaks=q.v, right=T, include.lowest=T, labels=PRS.labels)
phen$PRS.c <- relevel(PRS.c, ref = PRS.ref)


### Run logistic regression model adjusting for potential confounders
reg <- glm(t2d_status ~ PRS.c + age + sex + bmi + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = phen, family = binomial)
reg.s <- summary(reg)


