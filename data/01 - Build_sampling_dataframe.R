# Script to retrieve sampling data for each individual from the dataframe of sampling data for each cell
# It uses the names of the individuals, e.g. "C51i6" is SamplingCell #51, individual #6

rm(list=ls())

library(tidyverse)

# data.frame containing longitude, latitude, ecoregions, sampler etc for each samlping cell
load("Cell_sampling.RData")

# Diplodus sargus

# Read adaptive vcf (bcs it's smaller) and transform to genind
# (the genind format makes it easy to retrieve individual names
n_vcf <- read.vcfR("dip_adaptive_413.vcf.gz")
x <- vcfR2genind(n_vcf)

row.names(x@tab) # These are individual names

write.csv(row.names(x@tab), file="name_ind_Diplodus_adaptive.csv")

# in Excel, substr individual names to extract SamplingCcell, save into a csv file.

# Reload csv file
nomi <- read.csv2("name_ind_Diplodus_adaptive.csv")

# Build dataframe with individual names and sampling cell
Diplodus_sampling <- data.frame(ind_name = row.names(x@tab), SamplingCell = nomi$eccolo)

# Join sampling cell metadata to each individual
Diplodus_sampling <- left_join(Diplodus_sampling, cell_sampling, by = "SamplingCell")

Diplodus_sampling$SamplingCell <- factor(Diplodus_sampling$SamplingCell)

rm(n_vcf, x, nomi)

####


# Mullus surmuletus

# Read adaptive vcf (bcs it's smaller) and transform to genind
# (the genind format makes it easy to retrieve individual names
n_vcf <- read.vcfR("mul_adaptive_291.vcf.gz")
x <- vcfR2genind(n_vcf)

row.names(x@tab) # These are individual names

write.csv(row.names(x@tab), file="name_ind_Mullus_adaptive.csv")

# in Excel, substr individual names to extract SamplingCcell, save into a csv file.

# Reload csv file
nomi <- read.csv2("name_ind_Mullus_adaptive.csv")

# Build dataframe with individual names and sampling cell
Mullus_sampling <- data.frame(ind_name = row.names(x@tab), SamplingCell = nomi$eccolo)

# Join sampling cell metadata to each individual
Mullus_sampling <- left_join(Mullus_sampling, cell_sampling, by = "SamplingCell")

Mullus_sampling$SamplingCell <- factor(Mullus_sampling$SamplingCell)

rm(n_vcf, x, nomi)

###

save(cell_sampling, Diplodus_sampling, Mullus_sampling, file="sampling.RData")

