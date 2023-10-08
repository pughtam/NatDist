The anthropogenic imprint on temperate and boreal forest demography and carbon turnover

## Authors
Thomas A. M. Pugh, Lund University, Sweden, thomas.pugh@nateko.lu.se
Cornelius Senf, Technical University of Munich, corneliussenf@gmail.com

## Description
Scripts underlying analysis in:
Thomas A. M. Pugh, Rupert Seidl, Daijun Liu, Mats Lindeskog, Louise P. Chini, Cornelius Senf, The anthropogenic imprint on temperate and boreal forest demography and carbon turnover. Global Ecology and Biogeography.

Folder: disturbance_rates
- aggregate_hansen_forestcover.R, aggregate forest cover and disturbances from Hansen to definition of open/closed-canopy forests used by Pugh et al. 2019 Nature Geosciences
- prepare_hansen_disturbance_data.R, preparing Hansen disturbance data for non-spatial analysis, i.e., aggregating disturbance/forest pixels per year/landscape
- average_turover_rates.R, calibrate disturbance rate model for all forests
- average_turover_rates_closedforest.R, calibrate disturbance rate model for closed canopy forests
- cluster_assignment.R, map disturbance acitivity clusters to traits/climte using a multi-modal model

Folder: ESA_processing
- esa_broadleaf_frac_0p5deg.m, classify ESA landcover into broadleaf and needleleaf classes
- lpjg_foresttype_input.m, create an input file for LPJ-GUESS that specifies whether the forest type in each grid cell is broadleaf, needleleaf or mixed
- esa_lu_read.m, read and process the ESA landcover into general landcover types. The resulting landcover files is only used for the ocean mask in plots.

Folder: forest_mask
- hansen_forest_canopy_frac_calc.m, read in the forest cover from Hansen et al. (2013) and calculate the open and closed canopy forest cover, actual canopy cover and actual canopy area.

Folder: obs_biomes
- hengl_biome_read.m, Read Hengl et al. (2018) potential natural vegetation biome distribution
- hickler_biome_read.m, Read Haxeltine and Prentice (1996) potential natural vegetation biome distribution
- vegmap18.out, data file for Haxeltine and Prentice biomes

Folder: lpjg_processing
- trait_to_lpjg_mapping_tempbor.m, script to map trait values from species level to LPJ-GUESS PFTs (Figure S2, parameterisation for LPJ-GUESS)
- plot_global_distint_uncer_map.m, plot map of uncertainty in disturbance rotation periods, showing the absolute range of uncertainty divided by the best estimate (Figure S4).
- plot_global_cturn_diff_map.m, plot map of difference in carbon turnover from two different simulations and write processed data out to a text file (used for C turnover calculations in text).
- age_dist_plot.m, create age structure plots by region based on LPJ-GUESS output, including Figure S7 and Figure S8. Creates files used to draw Figure 4 using figure4.R.
- gfad_breakdown_region.m, split up the GFAD age class data into regions. Dependency for age_dist_plot.m
- readmasks_func.m, create the masks for low vegetation C mass, low forested area or non temperate/boreal biomes.
- hengl_lpjg_biome_comp.m, plot biomes from LPJ-GUESS in comparison to those from observation-based estimates (Figure S5).
- lpjg_biome_func.m, function to convert LPJ-GUESS LAI output into a biome classification.
- plot_global_distint_map.m, plot map of disturbance rotation periods and write processed data out to a text file (creates files used to draw Figure 3 using figure3.R).
- plot_pft_frac_map.m, script to make map of fraction of tree cover which is broadleaf for LPJ-GUESS and ESA data (Figure S6).
- distint_site_extract, extract disturbance return periods for specific sites to provide data for Table S4.
- compare_biomass_esa.m, extract biomass information from ESA CCI biomass for each of the 77 landscapes and compare to LPJ-GUESS results for the same locations (Figure S9).
- trait_space_map.m, for each grid cell in the LPJ-GUESS simulations, calculate whether it falls within the trait and climate space of the 77 landscapes (Figure S3).
- lpjg_succession_plots.m, make plot for succession tests (Figure S1).
- figure3.R, code for producing Figure 3 from files created by plot_global_distint_map.m
- figure4.R, code for producing Figure 4 from files created by age_dist_plot.m

Folder: lpjg_processing/helper_functions
- gfad_regions.m, script to take the GFAD regions and process them into something more accessible for plotting.
- reg_simplify.m, further simplify (i.e. aggregate) GFAD regions for plotting.
- global_grid_area.m, calculate the area of each grid cell.
- lpj_to_grid_func_centre.m, function to read in a LPJ-GUESS output file and reformat it into a multi-dimensional array suitable for making Matlab plots with.

Folder: lpjg_processing/netcdf_creation
- write_netcdf_lpjg.m, function to write raw LPJ-GUESS from gridded simulations to netcdfs for use by the rest of this processing chain.
- write_netcdf_lpjg_site.m, function to write raw LPJ-GUESS from site simulations to netcdfs for use by the rest of this processing chain.
- lpjg_to_netcdf.m, master function to read the LPJ-GUESS gridded data and call the routine to write it to netcdf
- lpjg_to_netcdf_site.m, master function to read the LPJ-GUESS site data and call the routine to write it to netcdf
Note: It is the netcdfs from these functions that are included in the associated data product, as they are self-documenting. Raw LPJ-GUESS text outputs used by the files in these outputs are not archived.

