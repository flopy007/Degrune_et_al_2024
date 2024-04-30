
# LOADING RAW DATA AND SUBSET CROPLAND SAMPLES --------------------------------

# load the otu tables  
otutab_cerco = read.table(
  "data/initiaData/otu_cerco.csv",
  h = T,
  sep = ";",
  row.names = 1
)
otutab_bact = read.table(
  "data/initialData/otu_bact.csv",
  h = T,
  sep = ";",
  row.names = 1
)
otutab_fung = read.table(
  "data/initialData/otu_fung.csv",
  h = T,
  sep = ";",
  row.names = 1
)
# load the taxonomy tables
taxo_cerco = read.table(
  "data/initialData/taxo_cerco.csv",
  h = T,
  sep = ";",
  row.names = 1
)
taxo_bact = read.table(
  "data/initialData/taxo_bact.csv",
  h = T,
  sep = ";",
  row.names = 1
)
taxo_fung = read.table(
  "data/initialData/taxo_fung.csv",
  h = T,
  sep = ";",
  row.names = 1
)
# load the meta data & management data
meta_data = read.table(
  "data/initialData/metaDATA.csv",
  h = T,
  sep = ";",
  row.names = 1
)
management_data = read.table(
  "data/initialData/management.csv",
  h = T,
  sep = ";",
  row.names = 1
)
# remove grassland samples and keep only cropland samples (n=156)
meta_data = meta_data[-which(meta_data$landuse_ID == "gr"), ]
otutab_cerco = otutab_cerco[, colnames(otutab_cerco) %in% rownames(meta_data)]
otutab_bact = otutab_bact[, colnames(otutab_bact) %in% rownames(meta_data)]
otutab_fung = otutab_fung[, colnames(otutab_fung) %in% rownames(meta_data)]
management_data = management_data[rownames(management_data) %in% rownames(meta_data), ]
meta_data$betadisp="1" # add a factor to compute distance from centroid later in the analysis

# export cleaned data as csv files
write.csv(otutab_cerco, file = "data/otutab_cercozoa.csv")
write.csv(otutab_bact, file = "data/otutab_bacteria.csv")
write.csv(otutab_fung, file = "data/otutab_fungi.csv")
write.csv(taxo_cerco, file = "data/taxonomy_cercozoa.csv")
write.csv(taxo_bact, file = "data/taxonomy_bacteria.csv")
write.csv(taxo_fung, file = "data/taxonomy_fungi.csv")
write.csv(meta_data, file = "data/meta_data.csv")
write.csv(management_data, file = "data/management_data.csv")


