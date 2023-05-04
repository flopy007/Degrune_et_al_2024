# This is code to replicate the analyses and figures from the paper Biodiversa
# Code developed by Florine Degrune and validated by Masahiro Ryo
# last update 16 March 2023


# setwd
setwd("[THE DIRECTORY WHERE YOU DOWNLOADED THE FOLDER]/Degrune_et_al_2023-main/")


rm(list=ls()) # to clear the global environment

# LOAD PACKAGES
pkgs = c("vegan", "data.table", "pastecs", "maps", "ggplot2", "tidyverse", "ggpubr", "lavaan", "mixOmics", "microbiome") # package names
inst = lapply(pkgs, library, character.only = TRUE) # load them
    #-------------------------------------------------------------------------------------------------------------------
    # NOTE: the package "mixOmics" was removed from CRAN, so one needs to install it by downloading the old tar.gz file 
    # from this website (https://cran.r-project.org/src/contrib/Archive/mixOmics/). 
    # Moreover, the package has some dependencies (e.g., igraph). These libraries also need to be installed.
    # For installing the package "microbiome", follow the instruction in https://microbiome.github.io/tutorials/Installation.html
    #-------------------------------------------------------------------------------------------------------------------


# SUBSET CROPLAND DATA --------------------------------------------------------
# data are subset to keep only cropland samples and remove grassland samples
# load the otu tables of cercozoa, bacteria and fungi 
otutab_cerco = read.table("data/initialData/otu_cerco.csv", h = T, sep = ";", row.names = 1)
otutab_bact = read.table("data/initialData/otu_bact.csv", h = T, sep = ";", row.names = 1)
otutab_fung = read.table("data/initialData/otu_fung.csv", h = T, sep = ";", row.names = 1)
# load the taxonomy tables of cercozoa, bacteria and fungi
taxo_cerco = read.table("data/initialData/taxo_cerco.csv", h = T, sep = ";", row.names = 1)
taxo_bact = read.table("data/initialData/taxo_bact.csv", h = T, sep = ";", row.names = 1)
taxo_fung = read.table("data/initialData/taxo_fung.csv", h = T, sep = ";", row.names = 1)
# load meta data
meta_data = read.table("data/initialData/metaDATA.csv", h = T, sep = ";", row.names = 1)
# select the variables (we removed bact16S and fungITS as they are not part of the hypothesis, and BD and MWD are removed too
# as there were no prior evidence on their impact on protists. We keep abiotic factors based on previous study)
meta_data = select(meta_data, -c(8:9, 18:19))
# remove grassland samples and keep only cropland samples (n=156)
meta_data = meta_data[-which(meta_data$landuse_ID == "gr"), ]
otutab_cerco = otutab_cerco[, colnames(otutab_cerco) %in% rownames(meta_data)]
otutab_bact = otutab_bact[, colnames(otutab_bact) %in% rownames(meta_data)]
otutab_fung = otutab_fung[, colnames(otutab_fung) %in% rownames(meta_data)]

# export data as csv files
write.csv(otutab_cerco, file = "data/otutab_cercozoa.csv")
write.csv(otutab_bact, file = "data/otutab_bacteria.csv")
write.csv(otutab_fung, file = "data/otutab_fungi.csv")
write.csv(taxo_cerco, file = "data/taxonomy_cercozoa.csv")
write.csv(taxo_bact, file = "data/taxonomy_bacteria.csv")
write.csv(taxo_fung, file = "data/taxonomy_fungi.csv")
write.csv(meta_data, file = "data/meta_data.csv")

# PREPARE DATA FOR ANALYSIS ------------------------------------
rm(list=ls()) # to clear the global environment
# 1. load cropland data ----
otutab_c = read.csv("data/otutab_cercozoa.csv", header = T, sep = ",", row.names = 1) # otu table cercozoa
otutab_b = read.csv("data/otutab_bacteria.csv", header = T, sep = ",", row.names = 1) # otu table bacteria
otutab_f = read.csv("data/otutab_fungi.csv", header = T, sep = ",", row.names = 1) # otu table fungi
taxonomy_c = read.csv("data/taxonomy_cercozoa.csv", header = T, sep = ",", row.names = 1) # taxonomy table protists
taxonomy_b = read.csv("data/taxonomy_bacteria.csv", header = T, sep = ",", row.names = 1) # taxonomy table bacteria
taxonomy_f = read.csv("data/taxonomy_fungi.csv", header = T, sep = ",", row.names = 1) # taxonomy table fungi
meta_data = read.csv("data/meta_data.csv", header = T, sep = ",", row.names = 1) # meta data

# 2. prepare microbial data ----
# cercozoa
otu = otutab_c
tax = taxonomy_c
source("scripts/dataPreparationCercozoa.R")
# bacteria
otu = otutab_b
tax = taxonomy_b
source("scripts/dataPreparationBacteria.R")
# fungi
otu = otutab_f
tax = taxonomy_f
source("scripts/dataPreparationFungi.R")
# combine data into a dataframe or list for later use
# data for RF and SEM analysis
data = cbind(meta_data,
             pco1_b,
             pco1_c,
             pco1_f,
             pco1_omnivore,
             pco1_bacterivore)
save(data, file="data/data.RData")      
# principal coordinates
pco_list = list(pco_b, pco_c, pco_f, pco_bacterivore, pco_omnivore)
names(pco_list) = c("pco_b", "pco_c", "pco_f", "pco_bacterivore", "pco_omnivore")
# bray curtis dissimilarity distance matrices
dist_list = list(dist_b, dist_c, dist_f, dist_bacterivore, dist_omnivore)
names(dist_list) = c("dist_b_dissim", "dist_c_dissim", "dist_f_dissim", "dist_bacteri_dissim", "dist_omni_dissim")
# taxonomy data
taxonomy_list = list(tax_b, tax_c, tax_f)
names(taxonomy_list) = c("tax_b", "tax_c", "tax_f")
# otu tables
otu_list = list(otu_b, otu_c, otu_f, otu_bacterivore, otu_omnivore)
names(otu_list) = c("otu_b", "otu_c", "otu_f", "otu_bacterivore", "otu_omnivore")
# save the lists in a RData file
save(pco_list,
     dist_list,
     taxonomy_list,
     otu_list,
     file = "data/microbial_lists.RData")

# 3. prepare taxonomy and functional data ----
rm(list=ls()) # to clear the global environment
# load the taxonomy data
load("data/microbial_lists.RData")
# cercozoa taxonomy and nutrition mode
# sum reads by FAMILY
df = taxonomy_list$tax_c[, c(4, 9:ncol(taxonomy_list$tax_c))]
df = aggregate(formula = . ~ family,
               data = df,
               FUN = sum)
df = data.frame(df[,-1], row.names = df[, 1]) # move first column to rowname
family_c = as.data.frame(t(df))
# compute relative abundance per sample
family_c_tot = decostand(family_c, method = "total")
# check if sum sample = 1
apply(family_c_tot, 1, sum)
# compute the percentage of reads for each family across all samples
sumseq = apply(family_c, 2, sum)
prop = (sumseq / sum(sumseq)) * 100
(family_c_prop = as.data.frame(t(rbind(sumseq, prop))))
# compute the sum of reads per nutrition mode
df = taxonomy_list$tax_c[, c(7, 9:ncol(taxonomy_list$tax_c))]
df = aggregate(formula = . ~ nutrition,
               data = df,
               FUN = sum)
df = data.frame(df[, -1], row.names = df[, 1]) # move first column to rowname
nutri_c = as.data.frame(t(df))
# compute relative abundance per sample
nutri_c_tot = decostand(nutri_c, method = "total")
# check if sum sample = 1
apply(nutri_c_tot, 1, sum)
# compute the percentage of reads for each family across all samples
sumseq = apply(nutri_c, 2, sum)
prop = (sumseq / sum(sumseq)) * 100
(nutri_c_prop = as.data.frame(t(rbind(sumseq, prop))))

# fungi taxonomy
# sum reads by PHYLUM
df = taxonomy_list$tax_f[, c(2, 8:ncol(taxonomy_list$tax_f))]
df = aggregate(formula = . ~ Phylum,
               data = df,
               FUN = sum)
df = data.frame(df[, -1], row.names = df[, 1]) # move first column to rowname
phylum_f = as.data.frame(t(df))
# compute relative abundance per sample
phylum_f_tot = decostand(phylum_f, method = "total")
# check if sum sample = 1
apply(phylum_f_tot, 1, sum)
# compute the percentage of reads for each family across all samples
sumseq = apply(phylum_f, 2, sum)
prop = (sumseq / sum(sumseq)) * 100
(phylum_f_prop = as.data.frame(t(rbind(sumseq, prop))))

# save taxonomy and functional data output into lists
taxa_countSample_list = list(family_c, phylum_f) # read count per sample
names(taxa_countSample_list) = c("family_c", "phylum_f")
taxa_relSample_list = list(family_c_tot, phylum_f_tot) # relative abundance per sample
names(taxa_relSample_list) = c("family_c_tot", "phylum_f_tot")
taxa_propOverall_list = list(family_c_prop, phylum_f_prop) # overall proportion of each taxa across all samples
names(taxa_propOverall_list) = c("family_c_prop", "phylum_f_prop")
nutrition_list = list(nutri_c, nutri_c_prop, nutri_c_tot)
names(nutrition_list) = c("nutri_c", "nutri_c_prop", "nutri_c_tot")

save(
  taxa_countSample_list,
  taxa_relSample_list,
  taxa_propOverall_list,
  nutrition_list,
  file = "data/taxonomy_lists.RData"
)

