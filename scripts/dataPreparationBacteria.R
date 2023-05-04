# combine taxonomy with otu table and keep only bacteria
tax1 = as.data.frame(cbind(tax, otu))
tax1 = tax1[which(tax1$kingdom == "Bacteria"),]
otu1 = as.data.frame(t(tax1[, 8:ncol(tax1)]))
# remove low abundant OTUs (<0.001%)
abd = (rowSums(tax1[, 8:ncol(tax1)]) / sum(tax1[, 8:ncol(tax1)])) * 100
ind1 = which(abd >= 0.001)
tax_b = tax1[ind1,]
otu_b = as.data.frame(t(tax_b[, 8:ncol(tax_b)]))
# standardization, compute distance matrix & ordinate data
otu_b_hel = decostand(otu_b, method = "hellinger")
dist_b = vegdist(otu_b_hel, method = "bray") # compute distance matrix using Bray-Curtis
pco_b = cmdscale(dist_b,
                 k = nrow(otu_b_hel) - 1,
                 eig = T,
                 add = T) # ordinate the data
# extract first principal component (pco1 which will be used in the random forest and SEM analysis)
pco1_b = pco_b$points[,1]
