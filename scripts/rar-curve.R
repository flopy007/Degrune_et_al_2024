

sample_names=row.names(otu)

rare <- lapply(out, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names

rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

data$sample=row.names(data)

rare <- map_dfr(data$sample, function(x){ #loop over samples
  z <- rare[rare$sample == x,] #subset rare according to sample 
  loc <- data$country[data$sample == x] #subset groupings according to sample, if more than one grouping repeat for all
  z <- data.frame(z, loc) #make a new data frame with the subsets
  return(z)
})

#### https://stackoverflow.com/questions/47234809/coloring-rarefaction-curve-lines-by-metadata-vegan-package-phyloseq-package




