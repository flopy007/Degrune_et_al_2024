# This is code to replicate the analyses and figures from the paper Biodiversa
# Code developed by Florine Degrune and validated by Masahiro Ryo
# last update 4th May 2023

# set the working directory
setwd("[THE DIRECTORY WHERE YOU DOWNLOADED THE FOLDER]/Degrune_et_al_2024-main/")

rm(list=ls()) # clear the global environment

# LOAD PACKAGES -----------------------------------------------------------------------
pkgs = c("vegan", "data.table", "pastecs", "maps", "ggplot2", "tidyverse", "ggpubr", "lavaan", "mixOmics", "RColorBrewer") # package names
inst = lapply(pkgs, require, character.only = TRUE) # load them

# ANALYSIS -----------------------------------------------------------------------
# 1. descriptive statistics OTUs and reads -------------
load("data/microbial_lists.RData") # LOAD DATA
# calculate descriptive statistics of cercozoan OTUs
(otu_sum=as.numeric(ncol(otu_list$otu_c)))
(otu_min=as.numeric(min(specnumber(otu_list$otu_c)))) 
(otu_max=as.numeric(max(specnumber(otu_list$otu_c))))
(otu_mean=as.numeric(mean(specnumber(otu_list$otu_c))))
(otu_median=as.numeric(median(specnumber(otu_list$otu_c))))
(otu_sd=as.numeric(sd(specnumber(otu_list$otu_c))))
sort(specnumber(otu_list$otu_c))
# calculate descriptive statistics of cercozoan reads
sort(rowSums(otu_list$otu_c))
(seq_sum=as.numeric(sum(otu_list$otu_c)))
(seq_min=as.numeric(min(rowSums(otu_list$otu_c))))
(seq_max=as.numeric(max(rowSums(otu_list$otu_c))))
(seq_mean=as.numeric(mean(rowSums(otu_list$otu_c))))
(seq_sd=as.numeric(sd(rowSums(otu_list$otu_c))))
(seq_median=as.numeric(median(rowSums(otu_list$otu_c))))

# calculate descriptive statistics of bacterial OTUs
(otu_sum=as.numeric(ncol(otu_list$otu_b)))
(otu_min=as.numeric(min(specnumber(otu_list$otu_b)))) 
(otu_max=as.numeric(max(specnumber(otu_list$otu_b))))
(otu_mean=as.numeric(mean(specnumber(otu_list$otu_b))))
(otu_median=as.numeric(median(specnumber(otu_list$otu_b))))
(otu_sd=as.numeric(sd(specnumber(otu_list$otu_b))))
sort(specnumber(otu_list$otu_b))
# calculate descriptive statistics of bacterial reads
sort(rowSums(otu_list$otu_b))
(seq_sum=as.numeric(sum(otu_list$otu_b)))
(seq_min=as.numeric(min(rowSums(otu_list$otu_b))))
(seq_max=as.numeric(max(rowSums(otu_list$otu_b))))
(seq_mean=as.numeric(mean(rowSums(otu_list$otu_b))))
(seq_sd=as.numeric(sd(rowSums(otu_list$otu_b))))
(seq_median=as.numeric(median(rowSums(otu_list$otu_b))))

# calculate descriptive statistics of fungal OTUs
(otu_sum=as.numeric(ncol(otu_list$otu_f)))
(otu_min=as.numeric(min(specnumber(otu_list$otu_f)))) 
(otu_max=as.numeric(max(specnumber(otu_list$otu_f))))
(otu_mean=as.numeric(mean(specnumber(otu_list$otu_f))))
(otu_median=as.numeric(median(specnumber(otu_list$otu_f))))
(otu_sd=as.numeric(sd(specnumber(otu_list$otu_f))))
sort(specnumber(otu_list$otu_f))
# calculate descriptive statistics of fungal reads
sort(rowSums(otu_list$otu_f))
(seq_sum=as.numeric(sum(otu_list$otu_f)))
(seq_min=as.numeric(min(rowSums(otu_list$otu_f))))
(seq_max=as.numeric(max(rowSums(otu_list$otu_f))))
(seq_mean=as.numeric(mean(rowSums(otu_list$otu_f))))
(seq_sd=as.numeric(sd(rowSums(otu_list$otu_f))))
(seq_median=as.numeric(median(rowSums(otu_list$otu_f))))

# calculate descriptive statistics of omnivorous OTUs
(otu_sum=as.numeric(ncol(otu_list$otu_omnivore)))
(otu_min=as.numeric(min(specnumber(otu_list$otu_omnivore)))) 
(otu_max=as.numeric(max(specnumber(otu_list$otu_omnivore))))
(otu_mean=as.numeric(mean(specnumber(otu_list$otu_omnivore))))
(otu_median=as.numeric(median(specnumber(otu_list$otu_omnivore))))
(otu_sd=as.numeric(sd(specnumber(otu_list$otu_omnivore))))
sort(specnumber(otu_list$otu_omnivore))
# calculate descriptive statistics of omnivorous reads
sort(rowSums(otu_list$otu_omnivore))
(seq_sum=as.numeric(sum(otu_list$otu_omnivore)))
(seq_min=as.numeric(min(rowSums(otu_list$otu_omnivore))))
(seq_max=as.numeric(max(rowSums(otu_list$otu_omnivore))))
(seq_mean=as.numeric(mean(rowSums(otu_list$otu_omnivore))))
(seq_sd=as.numeric(sd(rowSums(otu_list$otu_omnivore))))
(seq_median=as.numeric(median(rowSums(otu_list$otu_omnivore))))

# calculate descriptive statistics of bacterivorous cercozoa OTUs
(otu_sum=as.numeric(ncol(otu_list$otu_bacterivore)))
(otu_min=as.numeric(min(specnumber(otu_list$otu_bacterivore)))) 
(otu_max=as.numeric(max(specnumber(otu_list$otu_bacterivore))))
(otu_mean=as.numeric(mean(specnumber(otu_list$otu_bacterivore))))
(otu_median=as.numeric(median(specnumber(otu_list$otu_bacterivore))))
(otu_sd=as.numeric(sd(specnumber(otu_list$otu_bacterivore))))
sort(specnumber(otu_list$otu_bacterivore))
# calculate descriptive statistics of bacterivorous cercozoa reads
sort(rowSums(otu_list$otu_bacterivore))
(seq_sum=as.numeric(sum(otu_list$otu_bacterivore)))
(seq_min=as.numeric(min(rowSums(otu_list$otu_bacterivore))))
(seq_max=as.numeric(max(rowSums(otu_list$otu_bacterivore))))
(seq_mean=as.numeric(mean(rowSums(otu_list$otu_bacterivore))))
(seq_sd=as.numeric(sd(rowSums(otu_list$otu_bacterivore))))
(seq_median=as.numeric(median(rowSums(otu_list$otu_bacterivore))))

# 2. rarefaction curves and accumulation curve -------------
# rarefaction curves (sampling depth)
# apply rarecurve function
out= rarecurve(otu_list$otu_c, step=5, sample = 5000, xlab="Sample size", ylab="species", label=TRUE, col=data$country_ID)
save(out, file = "./output/rarefaction_curves/rar-curves-Cercozoa.RData")
out= rarecurve(otu_list$otu_b, step=5, sample = 5000, xlab="Sample size", ylab="species", label=TRUE, col=data$country_ID)
save(out, file = "./output/rarefaction_curves/rar-curves-Bacteria.RData")
out= rarecurve(otu_list$otu_f, step=5, sample = 5000, xlab="Sample size", ylab="species", label=TRUE, col=data$country_ID)
save(out, file = "./output/rarefaction_curves/rar-curves-Fungi.RData")
# plot rarefacion curves in a nice way
# cercozoa
load("./output/rarefaction_curves/rar-curves-Cercozoa.RData")
otu = otu_list$otu_c
source("scripts/rar-curve.R")
head(rare)
p.c=ggplot(data = rare, aes(x = max_raw, y = max_OTU, label=sample)) +
  geom_line(aes(x = raw.read, y = OTU, group = sample, color = loc))+
  scale_color_manual(values=c("#1F78B4" ,"#B2DF8A", "#E31A1C" ,"#FF7F00", "#6A3D9A")) +
  theme_bw()+theme(legend.position = c(0.8, 0.4))+
  labs(x="number of sequences - Cercozoa", y="number of OTUs - Cercozoa", col="Country") 
# bacteria
load("./output/rarefaction_curves/rar-curves-Bacteria.RData")
otu = otu_list$otu_b
source("scripts/rar-curve.R")
head(rare)
p.b=ggplot(data = rare, aes(x = max_raw, y = max_OTU, label=sample)) +
  geom_line(aes(x = raw.read, y = OTU, group = sample, color = loc))+
  scale_color_manual(values=c("#1F78B4" ,"#B2DF8A", "#E31A1C" ,"#FF7F00", "#6A3D9A")) +
  theme_bw()+theme(legend.position = c(0.8, 0.4))+
  labs(x="number of sequences - Bacteria", y="number of OTUs - Bacteria", col="Country") 
# fungi
load("./output/rarefaction_curves/rar-curves-Fungi.RData")
otu = otu_list$otu_f
source("scripts/rar-curve.R")
head(rare)
p.f=ggplot(data = rare, aes(x = max_raw, y = max_OTU, label=sample)) +
  geom_line(aes(x = raw.read, y = OTU, group = sample, color = loc))+
  scale_color_manual(values=c("#1F78B4" ,"#B2DF8A", "#E31A1C" ,"#FF7F00", "#6A3D9A")) +
  theme_bw()+theme(legend.position = c(0.8, 0.4))+
  labs(x="number of sequences - Fungi", y="number of OTUs - Fungi", col="Country") 

# plot rarefaction curves
svg("output/figure/rarefaction_curves.svg", width=12, height=4, family="sans")
ggarrange(p.c, p.b, p.f, ncol = 3, nrow = 1)
dev.off()

# accumulation curve (sampling effort)
svg("output/figure/accumulation_curve.svg", width=8, height=3, family="sans")
par(mfrow = c(1, 3))
plot(specaccum(otu_list$otu_c), ylab="number of cercozoan OTUs", xlab="number of samples")
plot(specaccum(otu_list$otu_b), ylab="number of bacterial OTUs", xlab="number of samples")
plot(specaccum(otu_list$otu_f), ylab="number of fungal OTUs", xlab="number of samples")
dev.off()


# 3. random forest analysis ----
rm(list=ls()) # CLEAR THE GLOBAL ENVIRONMENT
load("data/data.RData") # LOAD DATA
# create directories
dir.create(path = "output/modeling_rf/RF_cercozoa_pco1")
dir.create(path = "output/modeling_rf/RF_cercozoa_pco2")
dir.create(path = "output/modeling_rf/RF_omnivore_pco1")
dir.create(path = "output/modeling_rf/RF_omnivore_pco2")
dir.create(path = "output/modeling_rf/RF_bacterivore_pco1")
dir.create(path = "output/modeling_rf/RF_bacterivore_pco2")
# the computing time is long. Consider a few hours on a regular machine.
source("scripts/RF_cercozoa_pco1.R") # random forest to predict Cercozoa PCo1
source("scripts/RF_cercozoa_pco2.R") # random forest to predict Cercozoa PCo2
source("scripts/RF_omnivore_pco2.R") # random forest to predict omnivore PCo2
source("scripts/RF_omnivore_pco1.R") # random forest to predict omnivore PCo1
source("scripts/RF_bacterivore_pco1.R") # random forest to predict bacterivore pco1
source("scripts/RF_bacterivore_pco2.R") # random forest to predict bacterivore pco2

# the output can be found in the different folders in output/modeling_rf
# proceed manually to format the output file
# the formated file is available here: output/modeling_rf/rf_varImp_output.csv
# load the formated file
rf_output = read.csv("output/modeling_rf/rf_varImp_output.csv", h = T, sep = ",")
colnames(rf_output)[1] = "Name"
# reshape table
table_rf = reshape2::melt(rf_output)
table_rf_1 = table_rf[which(table_rf$Name == "Cercozoa_pco1" |  table_rf$Name == "Cercozoa_pco2"),]
table_rf_2 = table_rf[which(table_rf$Name == "Omnivore_pco1" |  table_rf$Name == "Bacterivore_pco1"),]
# plot 1
svg("./output/figure/barplot_RF_cercozoa_pco1_pco2.svg", 
    width=4, height=4, family="sans")
(p<-ggplot(data=table_rf_1, aes(x=value, y=variable, fill = Name)) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_fill_manual(values=c("#B15928", "#FDBF6F"))+
    theme_minimal() +
    theme(legend.position="top"))
dev.off()
# plot 2
svg("./output/figure/barplot_RF_omnivore_bacterivore_pco1.svg", 
    width=4, height=4, family="sans")
(p<-ggplot(data=table_rf_2, aes(x=value, y=variable, fill = Name)) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_fill_manual(values=c("#B2DF8A" ,"#33A02C"))+
    theme_minimal() +
    theme(legend.position="top"))
dev.off()

# open the svg file with inkscape and make esthetic modifications

# 4. structural equation modeling (SEM) ----
rm(list=ls()) # CLEAR THE GLOBAL ENVIRONMENT
load("data/data.RData") # LOAD DATA
# store paths
output_file = "./output/modeling_sem/output_SEM.csv"
# subset data for buliding the SEM (abiotic factors that are significant in the abiotic-biotic RF model)
colnames(data)
data_subset = data[, -c(1:5)] # remove metadata
colnames(data_subset)
# scale data
data_scale=as.data.frame(scale(data_subset))
colnames(data_scale)

# Model building
# Model 1
model <- '
pco1_c ~  MAP 
          + MAT
          + Lat
          + Long
          + SM
          + Clay 
          + Ctot_mg
          + Ntot_mg
          + Ptot_mg
          + pH
          + pco1_f
pco1_f ~ MAP 
          + MAT
          + Lat
          + Long
          + SM
          + Clay 
          + Ctot_mg
          + Ntot_mg
          + Ptot_mg
          + pH
'
set.seed(3463)
fit = lavaan(model, data=data_scale, auto.var=TRUE)
summary(fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
lavaan::anova(fit)

# Model 2 (remove the non-significant link from model 1)
model <- '
pco1_c ~  MAP 
         # + MAT
          + Lat
          #+ Long
          + SM
          + Clay 
          #+ Ctot_mg
          #+ Ntot_mg
          + Ptot_mg
          + pH
          + pco1_f
pco1_f ~ MAP 
         # + MAT
          + Lat
          #+ Long
          #+ SM
          + Clay 
          + Ctot_mg
          + Ntot_mg
          + Ptot_mg
          #+ pH
'
set.seed(3463)
fit = lavaan(model, data=data_scale, auto.var=TRUE)
summary(fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
lavaan::anova(fit)

# Model 3 (remove the non-significant link from model 2)
model <- '
pco1_c ~  MAP 
         # + MAT
          + Lat
          #+ Long
          #+ SM
          + Clay 
          #+ Ctot_mg
          #+ Ntot_mg
          + Ptot_mg
          + pH
          + pco1_f
pco1_f ~ MAP 
         # + MAT
          + Lat
          #+ Long
          #+ SM
          + Clay 
          + Ctot_mg
          + Ntot_mg
          + Ptot_mg
          #+ pH
'
set.seed(3463)
fit = lavaan(model, data=data_scale, auto.var=TRUE)
summary(fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
lavaan::anova(fit)

# generate results
results = parameterEstimates(fit, standardized = T, rsquare = TRUE)
# export results
write.csv(results, output_file)

# 5. distance decay analysis --------------
rm(list=ls()) # CLEAR THE GLOBAL ENVIRONMENT
# LOAD DATA
load("data/microbial_lists.RData")
load("data/data.RData")
# extract geographic distances
geo=data[, c(12:13)] # extract latitude and longitude
geo$x.km=geo$Lat*111.000 # transform into geographic distance (km)
# compute a distance matrix based on euclidean distances
dist_geo=vegdist(geo$x.km, method="euclidean") # geographic dissimilarity
# extract the bray curtis dissimilarity of Cercozoa and fungi
dist_f = dist_list$dist_f_dissim # dissimilarity in community composition
dist_c = dist_list$dist_c_dissim # dissimilarity in community composition
## run linear regression
(dd.lm.f.geo = lm(dist_f ~ dist_geo)) # intercept = 6.327e-01, slope = 8.362e-05  
(dd.lm.c.geo = lm(dist_c ~ dist_geo)) # intercept = 4.532e-01 , slope = 5.053e-05
## summary & regression line
summary(dd.lm.f.geo) #### Adjusted R-squared: 0.3058  ; p-value: < 2.2e-16
summary(dd.lm.c.geo) #### Adjusted R-squared:  0.2226  ; pval< 2.2e-16
## plot
svg("./output/figure/dd_cercozoa.svg", width=5, height=5, family="sans")
plot(
  dist_c ~ dist_geo,
  col = "grey",
  type = "p",
  pch = 19,
  cex = 0.5,
  xlab = "geographic distance (km)",
  ylab = "Cercozoan community dissimilarity"
)
abline(4.532e-01, 5.053e-05, col="black", lwd = 2)
dev.off()

svg("./output/figure/dd_fungi.svg", width=5, height=5, family="sans")
plot(
  dist_f ~ dist_geo,
  col = "grey",
  type = "p",
  pch = 19,
  cex = 0.5,
  xlab = "geographic distance (km)",
  ylab = "Fungal community dissimilarity"
)
abline(6.327e-01, 8.362e-05 , col="black", lwd = 2)
dev.off()

# 6. linear multivariate regression (spls) ----
rm(list=ls()) # CLEAR THE GLOBAL ENVIRONMENT
# LOAD DATA
load(file="data/taxonomy_lists.RData")
# define x (predictor = fungi) and y (variable to predict = cercozoa)
# Y = Cercozoa at the family level, X = fungi at the phylum level
Y = taxa_relSample_list$family_c_tot
X_f = taxa_relSample_list$phylum_f_tot
# subset only omnivores for which feeding activity on fungi was observed or suspected (see Dumack et al., 2019) 
# and exclude Polymyxa_lineage and Spongospora_lineaage as they are known as plant parasites mostly
names_c = c("Assulinidae","Cercomonadidae","Fiscullidae","Plasmodiophoridae","Trinematidae",
            "Viridiraptoridae","Euglyphidae" , "Rhogostomidae","Sphenoderiidae",
            "Tremulidae")
Y1 = Y[, names_c]

# initial sPLS model
spls_initial = spls(X = X_f, Y = Y1, ncomp = ncol(X_f), mode = 'regression')
# tuning sPLS (selection the optimal number of component)
# repeated CV tuning of component count
perf_spls <- perf(spls_initial, validation = 'Mfold',
                    folds = 10, nrepeat = 5) 
abline(h = 0.0975)
perf_spls$measures$Q2.total$values
perf_spls$measures$R2$values
## we keep the two first components
# here, we want to keep all X variables in the model, we don't do variable selection
# final model
set.seed(3958)
spls_final = spls(X_f, Y1, ncomp = 2, 
                    mode = "regression") # explanatory approach being used, hence use regression mode
# plot Heat map correlation
colfunc = colorRampPalette(c("#B15928", "white", "#33A02C"))
svg("output/figure/spls_analysis.svg", width = 6, height = 6)
spls_final_plot = cim(spls_final, comp = 1:2, xlab = "", ylab = "", margins = c(4, 8), color=colfunc(20))
dev.off()
# output table with correlation values
matcor_output = as.data.frame(spls_final_plot$mat.cor)
matcor_output = setDT(matcor_output, keep.rownames = TRUE)[]
names(matcor_output)[1] = "clades"
row.names(matcor_output) = matcor_output$clades
matcor_output2 = matcor_output[, -1]
matcor_output2 = round(matcor_output2, digits = 3)
row.names(matcor_output2) = row.names(matcor_output)
write.table(matcor_output2, "output/modeling_spls/spls_analysis.csv", sep = ",") # export the table
# 7. ordination PCoA  ---------------------------------------------
rm(list=ls()) # CLEAR THE GLOBAL ENVIRONMENT
# LOAD DATA
load("data/data.RData") # load data
load("data/microbial_lists.RData")
## principal coordinate anlaysis Cercozoa
# get principal coordinates
cp = pco_list$pco_c$eig / sum(pco_list$pco_c$eig) * 100
cp = round(cp[1:2], digits = 2)### 11.44 / 9.04
# plot the ordination
svg("output/figure/pcoa_cercozoa.svg", width = 5, height = 5, family = "sans")
fig = ordiplot(
  pco_list$pco_c,
  type = "n",
  xlab = 'PCo1 - 11.44%',
  ylab = 'PCo2 - 9.04%',
  cex.axis = 1.0,
  cex.lab = 1.0,
  choices = c(1:2)
)
abline(fig, h = 0, v = 0, col = "honeydew3", lty = 3)
col = c("#1F78B4" ,"#B2DF8A", "#E31A1C" ,"#FF7F00", "#6A3D9A")[as.numeric(as.factor(data$country_ID))]
points(fig, "sites", pch = 16, col = col, cex = 0.8, lwd = 1)
# text(fig, "sites", labels=data$country_ID, col="gray47", cex=0.5, pos=1, family="sans")
# add ordihull
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "SWE")], col = "#FF7F00")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "GER")], col = "#B2DF8A")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "SWI")], col = "#6A3D9A")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "FRA")], col = "#1F78B4")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "SPA")], col = "#E31A1C")

legend("topright", legend=c("Sweden", "Germany", "Switzerland", "France", "Spain"), 
       bty="n", cex=0.8, horiz = F, 
       text.col=c("#FF7F00", "#B2DF8A", "#6A3D9A", "#1F78B4", "#E31A1C"),
       xpd=NA, inset = c(0.0004, 0.00001), xjust=0)
dev.off()

## principal coordinate anlaysis Fungi
## get principal coordinate
cp = pco_list$pco_f$eig / sum(pco_list$pco_f$eig) * 100
cp = round(cp[1:2], digits = 2)### 13.5 / 9.7
# plot the ordination
svg("output/figure/pcoa_fungi.svg", width = 5, height = 5, family = "sans")
fig = ordiplot(
  pco_list$pco_f,
  type = "n",
  xlab = 'PCo1 - 13.5%',
  ylab = 'PCo2 - 9.7%',
  cex.axis = 1.0,
  cex.lab = 1.0,
  choices = c(1:2)
)
abline(fig, h = 0, v = 0, col = "honeydew3", lty = 3)
col = c("#1F78B4" ,"#B2DF8A", "#E31A1C" ,"#FF7F00", "#6A3D9A")[as.numeric(as.factor(data$country_ID))]
points(fig, "sites", pch = 16, col = col, cex = 0.8, lwd = 1)
#text(fig, "sites", labels=substr(rownames(otu_c), 4, 6), col="gray47", cex=0.5, pos=1, family="sans")
# add ordihull
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "SWE")], col = "#FF7F00")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "GER")], col = "#B2DF8A")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "SWI")], col = "#6A3D9A")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "FRA")], col = "#1F78B4")
ordihull(fig, groups=data$country_ID, show.groups = data$country_ID[which(data$country_ID == "SPA")], col = "#E31A1C")

legend("topright", legend=c("Sweden", "Germany", "Switzerland", "France", "Spain"), 
       bty="n", cex=0.8, horiz = F, 
       text.col=c("#FF7F00", "#B2DF8A", "#6A3D9A", "#1F78B4", "#E31A1C"),
       xpd=NA, inset = c(0.0004, 0.00001), xjust=0)
dev.off()
