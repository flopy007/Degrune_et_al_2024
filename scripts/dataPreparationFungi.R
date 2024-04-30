# mapping taxonomy file with otutab
tax1 <- tax[rownames(tax) %in% rownames(otu), ]
# combine taxonomy with otu table and keep only fungi
tax2 = as.data.frame(cbind(tax1, otu))
tax2 = tax2[which(tax2$Kingdom == "Fungi"), ]
otu1 = as.data.frame(t(tax2[, 8:ncol(tax2)]))
# remove low abundant OTUs (<0.001% - see Degrune et al., 2019)
abd = (rowSums(tax2[, 8:ncol(tax2)]) / sum(tax2[, 8:ncol(tax2)])) * 100
ind1 = which(abd >= 0.001)
tax_f = tax2[ind1,]
otu_f = as.data.frame(t(tax_f[, 8:ncol(tax_f)]))
# standardization, compute distance matrix & ordinate data
otu_f_hel = decostand(otu_f, method = "hellinger")
dist_f = vegdist(otu_f_hel, method = "bray")
pco_f = cmdscale(dist_f,
                 k = nrow(otu_f_hel) - 1,
                 eig = T,
                 add = T)
# extract first principal component (pco1 which will be used in the random forest and SEM analysis)
pco1_f = pco_f$points[,1]
# plot principal coordinate
plot(pco_b$eig)
# we keep the two first components that summarize most of the variance
# extract the two first principal component (pco1 which will be used in the random forest and SEM analysis)
pco1_f = pco_f$points[,1]
pco2_f = pco_f$points[,2]

