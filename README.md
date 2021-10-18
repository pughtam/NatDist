The anthropogenic imprint on temperate and boreal forest demography and carbon turnover

## Authors
Thomas A. M. Pugh, Lund University, Sweden, thomas.pugh@nateko.lu.se
Cornelius Senf, Technical University of Munich, corneliussenf@gmail.com


## Description
Scripts underlying analysis in:
Thomas A. M. Pugh, Rupert Seidl, Daijun Liu, Mats Lindeskog, Louise P. Chini, Cornelius Senf, The anthropogenic imprint on temperate and boreal forest demography and carbon turnover.

Folder: disturbance_rates
- aggregate_hansen_forestcover.R, aggregate forest cover and disturbances from Hansen to definition of open/closed-canopy forests used by Pugh et al. 2019 Nature Geosciences
- prepare_hansen_disturbance_data.R, preparing Hansen disturbance data for non-spatial analysis, i.e., aggregating disturbance/forest pixels per year/landscape
- average_turover_rates.R, calibrate disturbance rate model for all forests
- average_turover_rates_closedforest.R, calibrate disturbance rate model for closed canopy forests
- cluster_assignment.R, map disturbance acitivity clusters to traits/climte using a multi-modal model

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
- trait_to_lpjg_mapping_tempbor.m, script to map trait values from species level to LPJ-GUESS PFTs (Figure S4, parameterisation for LPJ-GUESS)
- plot_global_distint_diff_map.m, plot map of difference between disturbance rotation periods from two different simulations and write processed data out to a text file (feeds into Figure 2).
- plot_global_distint_uncer_map.m, plot map of uncertainty in disturbance rotation periods, showing the absolute range of uncertainty divided by the best estimate (Figure S2, feeds into Figure 2).
- plot_global_cturn_diff_map.m, plot map of difference in carbon turnover from two different simulations and write processed data out to a text file (used for C turnover calculations in text).
- age_dist_plot.m, create age structure plots by region based on LPJ-GUESS output, including Figure S5 and Figure S6.
- gfad_breakdown_region.m, split up the GFAD age class data into regions. Dependency for age_dist_plot.m
- readmasks_func.m, create the masks for low vegetation C mass, low forested area or non temperate/boreal biomes.
- hengl_lpjg_biome_comp.m, plot biomes from LPJ-GUESS in comparison to those from observation-based estimates (Figure S1).
- lpjg_biome_func.m, function to convert LPJ-GUESS LAI output into a biome classification.
- plot_global_distint_map.m, plot map of disturbance rotation periods and write processed data out to a text file.
- plot_pft_frac_map.m, script to make map of fraction of tree cover which is broadleaf for LPJ-GUESS and ESA data (Figure S2).
- figure2.R, code for producing Figure 2
- figure3.R, code for producing Figure 3

Folder: lpjg_processing/helper_functions
- gfad_regions.m, script to take the GFAD regions and process them into something more accessible for plotting.
- reg_simplify.m, further simplify (i.e. aggregate) GFAD regions for plotting.
- global_grid_area.m, calculate the area of each grid cell.
- lpj_to_grid_func_centre.m, function to read in a LPJ-GUESS output file and reformat it into a multi-dimensional array suitable for making Matlab plots with.


