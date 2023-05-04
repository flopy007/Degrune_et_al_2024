# mapping taxonomy file with otutab
otu1 <- otu[rownames(otu) %in% rownames(tax), ]
# combine taxonomy with otu table
tax1 = as.data.frame(cbind(tax, otu1))
otu2 = as.data.frame(t(tax1[, 9:ncol(tax1)])) # samples as rows and otus as columns
# remove low abundant OTUs (<0.01% - see Fiore Donno et al., 2019)
abd = (rowSums(tax1[, 9:ncol(tax1)]) / sum(tax1[, 9:ncol(tax1)])) * 100
ind1 = which(abd >= 0.01)
tax_c = tax1[ind1,]
otu_c = as.data.frame(t(tax_c[, 9:ncol(tax_c)]))
# standardization, compute distance matrix & ordinate data - for all OTUs
otu_c_hel = decostand(otu_c, method = "hellinger")
dist_c = vegdist(otu_c_hel, method = "bray") # compute distance matrix using Bray-Curtis
pco_c = cmdscale(dist_c,
                 k = nrow(otu_c) - 1,
                 eig = T,
                 add = T) # ordinate the data
# standardization, compute distance matrix & ordinate data - for bacterivores OTUs
otu_bacterivore = as.data.frame(t(tax_c[which(tax_c$nutrition == "bacterivore"), 9:ncol(tax_c)]))
otu_bacterivore_hel = decostand(otu_bacterivore, method = "hellinger")
dist_bacterivore = vegdist(otu_bacterivore_hel, method = "bray") # compute distance matrix using Bray-Curtis
pco_bacterivore = cmdscale(dist_bacterivore,
                           k = nrow(otu_bacterivore) - 1,
                           eig = T,
                           add = T) # ordinate the data
# standardization, compute distance matrix & ordinate data - for omnivores OTUs
otu_omnivore = as.data.frame(t(tax_c[which(tax_c$nutrition == "omnivore"), 9:ncol(tax_c)]))
otu_omnivore_hel = decostand(otu_omnivore, method = "hellinger")
dist_omnivore = vegdist(otu_omnivore_hel, method = "bray") # compute distance matrix using Bray-Curtis
pco_omnivore = cmdscale(dist_omnivore,
                        k = nrow(otu_omnivore) - 1,
                        eig = T,
                        add = T) # ordinate the data
# extract first principal component (pco1 which will be used in the random forest and SEM analysis)
pco1_c = pco_c$points[,1]
pco1_bacterivore = pco_bacterivore$points[,1]
pco1_omnivore = pco_omnivore$points[,1]
