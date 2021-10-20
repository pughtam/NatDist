%Write netcdf files for all uptake data underlying manuscript figures for data deposition.
%The basis is tslice files from LPJ-GUESS which give means over the years 2001-2014.
%
%Dependencies: write_netcdf_lpjg.m
%
%T. Pugh
%20.10.21

input_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results';
output_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';

sims={'best_est_adjparam_latosa4000','low_2se_adjparam_latosa4000','high_2se_adjparam_latosa4000',...
'best_est_adjparam_latosa4000_luh2','best_est_adjparam_latosa4000_luh2low','best_est_adjparam_latosa4000_luh2high',...
'best_est_adjparam_latosa4000_closedcan'};

dtype={'nat','nat','nat',...
    'anthro','anthro','anthro'...
    'natcc'};
dreal={'standard','low','high',...
    'standard','low','high',...
    'standard'};
nsims=length(sims);

PFTs={'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G'};

if (length(dtype)~=nsims || length(dreal)~=nsims)
    error('Arrays describing simulations are not consistent in length')
end

%First write files for the simulations with age forcing for the full period
for ss=1:nsims
    
    currdir=[input_dir,'/',sims{ss},'/'];
    
    cd(currdir)
    
    %Get cpool data
    cpool=squeeze(lpj_to_grid_func_centre('cpool_2001_2014',1,0));
    cveg=cpool(:,:,1);
    csoil=cpool(:,:,3);
    clitter=cpool(:,:,2);
    clear cpool
    
    %Get cflux data
    cflux=squeeze(lpj_to_grid_func_centre('cflux_2001_2014',1,0));
    npp=cflux(:,:,1);
    clear cflux
    
    %Get LAI data
    lai=squeeze(lpj_to_grid_func_centre('lai_2001_2014',1,0));
    lai=lai(:,:,1:length(PFTs));
    
    %Write each dataset to netcdf file
    
    %Cveg
    variable='Cveg';
    variable_longname='Carbon stock in live vegetation';
    units='kg C m-2';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';
    
    write_netcdf_lpjg(cveg,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1);
    
    %Csoil
    variable='Csoil';
    variable_longname='Carbon stock in soil';
    units='kg C m-2';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';
    
    write_netcdf_lpjg(csoil,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1);
    
    %Clitter
    variable='Clitter';
    variable_longname='Carbon stock in litter';
    units='kg C m-2';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';
    
    write_netcdf_lpjg(clitter,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1);
    
    %NPP
    variable='NPP';
    variable_longname='Net primary productivity';
    units='kg C m-2 yr-1';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';
    
    write_netcdf_lpjg(npp,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1);
    
    %LAI
    variable='LAI';
    variable_longname='Leaf area index';
    units='m2 (leaf) m-2 (ground)';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';
    
    write_netcdf_lpjg(lai,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1,PFTs);
    
end
clear ss
