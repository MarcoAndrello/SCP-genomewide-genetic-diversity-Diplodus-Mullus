library(vcfR)
library(adegenet)

# Diplodus
# a_vcf <- read.vcfR("dip_neutral_7655.vcf.gz")
# b_vcf <- read.vcfR("dip_adaptive_413.vcf.gz")
# Mullus
a_vcf <- read.vcfR("mul_neutral_2462.vcf.gz")
b_vcf <- read.vcfR("mul_adaptive_291.vcf.gz")
a <- vcfR2genind(a_vcf)
b <- vcfR2genind(b_vcf)

# Same individual names?
all(row.names(a@tab) == row.names(b@tab))

# Get genotypes and individual names from both genind objects
genotypes_a <- a$tab
genotypes_b <- b$tab

# Combine the genotypes for the same individuals
merged_genotypes <- cbind(genotypes_a, genotypes_b)

# Create a new genind object with the merged genotypes
# Diplodus_8068 <- new("genind", tab = merged_genotypes)
Mullus_2753 <- new("genind", tab = merged_genotypes)

# save(Diplodus_8068,file="Diplodus_8068.RData")
save(Mullus_2753,file="Mullus_2753.RData")

