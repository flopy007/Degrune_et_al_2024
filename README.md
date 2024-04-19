# Degrune_et_al_2024

This folder contains the data and R codes to reproduce the data analysis and figures of the paper Degrune et al., 2023.

The **scripts** folder contain all R scripts needed to reproduce the analysis. You should start with **01_dataPreparation.R** and continue with **02_dataAnalysis**.
The unnumbered scripts are imbeded in the two R scripts and will be automatically run with the command source().

The **data** folder contains all files needed to run the analysis in R studio. The files in the data folder are of two types: the **initial** data (in the **initialData** folder) and data resulting from computing steps during the analysis (in the **data** folder). Therefore, if you want to start from scratch, I recommend to only keep the files in the **initialData** folder. The initial data can be found in the original publication *Garland et al., 2021* (see reference below).

The *output* folder contains the figures produced during the analysis (*output/figure*), the output of the random forest modeling (*output/modeling_rf*), the output of the structural equation modeling (*output/modeling_sem*), the output of the spls analysis (*output/modeling_spls*), and the output of the rarefaction curves (*output/rarefaction_curves*). The supplementary material of the mansucript is also provided in the output folder (*output/supplementary_material*)

For the figures, I usually save them as *.svg* and use the free and open source image editor **Inkscape** (version used: 1.0.2) (https://inkscape.org/fr/) to manipulate the image and save it as *.tif* or another format. The *.png* figures are then combined into a global *.xcf* figure using the free and open source image editor **GIMP** (version 2.10) (https://www.gimp.org). The resulted figure is then saved as *.png*. The main figures fig1 and fig2 of the manuscript are provided as *.xcf* files and *.tif* files in the folder *output/figure*. 

For the SEM analysis, the figure was manually created with Inkscape. The values in the figure (pvalue and the coefficient std.lv) are from the output_SEM.csv file located in the output folder > modeling_sem.

For more information on the way data were collected, please refer to Garland et al., 2021. Crop Cover Is More Important than Rotational Diversity for Soil Multifunctionality and Cereal Yields in European Cropping Systems. Nature Food 2 (1): 28–37. https://doi.org/10.1038/s43016-020-00210-8.

For more informations, contact me at florine.degrune@cirad.fr
