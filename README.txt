The anthropogenic imprint on temperate and boreal forest demography and carbon turnover

## Authors
Thomas A. M. Pugh, Lund University, Sweden, thomas.pugh@nateko.lu.se
Cornelius Senf, Technical University of Munich, corneliussenf@gmail.com


## Description
Scripts underlying analysis in:
Thomas A. M. Pugh, Rupert Seidl, Daijun Liu, Mats Lindeskog, Cornelius Senf, The anthropogenic imprint on temperate and boreal forest demography and carbon turnover.

Folder: disturbance_rates
- aggregate_hansen_forestcover.R, aggregate forest cover and disturbances from Hansen to definition of open/closed-canopy forests used by Pugh et al. 2019 Nature Geosciences
- prepare_hansen_disturbance_data.R, preparing Hansen disturbance data for non-spatial analysis, i.e., aggregating disturbance/forest pixels per year/landscape
- average_turover_rates.R, calibrate disturbance rate model for all forests
- average_turover_rates_closedforest.R, calibrate disturbance rate model for closed canopy forests
- cluster_assignment.R, map disturbance acitivity clusters to traits/climte using a multi-modal model
- figure2.R, code for producing Figure 2
- figure3.R, code for producing Figure 3

Folder: ESA_processing
- esa_broadleaf_frac_0p5deg.m, classify ESA landcover into broadleaf and needleleaf classes
- lpjg_foresttype_input.m, Create an input file for LPJ-GUESS that specifies whether the forest type in each grid cell is broadleaf, needleleaf or mixed

Folder: forest_mask
- hansen_forest_canopy_frac_calc.m, read in the forest cover from Hansen et al. (2013) and calculate the open and closed canopy forest cover, actual canopy cover and actual canopy area.

Folder: obs_biomes
- hengl_biome_read.m, Read Hengl et al. (2018) potential natural vegetation biome distribution
- hickler_biome_read.m, Read Haxeltine and Prentice (1996) potential natural vegetation biome distribution
- vegmap18.out, data file for Haxeltine and Prentice biomes

Folder: lpjg_processing

