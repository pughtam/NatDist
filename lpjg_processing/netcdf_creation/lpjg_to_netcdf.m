% Write netcdf files for all data underlying manuscript figures for data deposition.
% All other calculations with LPJ-GUESS data should be based on these netcdf files
% The basis is tslice files from LPJ-GUESS which give means over the years 2001-2014.
%
% Dependencies: write_netcdf_lpjg.m
%
% T. Pugh
% 20.10.21

addpath('./')

input_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results';
output_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';

luh2_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/NatDist_working_Jul23/lpjg_processing/data/';
luh2_file='lu_1700_2015_luh2_aggregate_sum2x2_midpoint_urban_orig_v20_2001_2014';

sims={'best_est_adjparam_IBS150BNE300_p100','low_2se_adjparam_IBS150BNE300','high_2se_adjparam_IBS150BNE300',...
'best_est_adjparam_IBS150BNE300_luh2','best_est_adjparam_IBS150BNE300_luh2low','best_est_adjparam_IBS150BNE300_luh2high',...
'best_est_adjparam_IBS150BNE300_closedcan'};

dtype={'nat','nat','nat',...
    'anthro','anthro','anthro'...
    'natcc'};
dreal={'standard','low','high',...
    'standard','low','high',...
    'standard'};
nsims=length(sims);

PFTs={'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G'};
ageclasses={'1-10','11-20','21-30','31-40','41-50','51-60','61-70','71-80','81-90','91-100','101-110','111-120','121-130','131-140','>140'};

if (length(dtype)~=nsims || length(dreal)~=nsims)
    error('Arrays describing simulations are not consistent in length')
end

%First write files for the simulations with age forcing for the full period
for ss=1:nsims
    
    currdir=[input_dir,'/',sims{ss},'/'];
    
    cd(currdir)
    
    if strcmp(dtype{ss},'anthro')
        cpool_file='cpool_natural_2001_2014';
        cflux_file='cflux_natural_2001_2014';
    else
        cpool_file='cpool_2001_2014';
        cflux_file='cflux_2001_2014';
    end
    
    %Get cpool data
    cpool=squeeze(lpj_to_grid_func_centre(cpool_file,1,0));
    cveg=cpool(:,:,1);
    csoil=cpool(:,:,3);
    clitter=cpool(:,:,2);
    clear cpool
    
    %Get cflux data
    cflux=squeeze(lpj_to_grid_func_centre(cflux_file,1,0));
    npp=-cflux(:,:,1);
    clear cflux

    %Get GPP data
    if ss==1
        agpp=squeeze(lpj_to_grid_func_centre('agpp_2001_2014',1,0));
        gpp=agpp(:,:,13);
        fprintf('Note: If changing code, need to take care to check correct column is used for GPP\n')
        clear agpp
    elseif ss==4
        agpp=squeeze(lpj_to_grid_func_centre('agpp_2001_2014',1,0));
        gpp=agpp(:,:,26);
        fprintf('Note: If changing code, need to take care to check correct column is used for GPP\n')
        clear agpp
        % Correct by LUH2 fraction here to keep consistent with other output files
        luh2_in=squeeze(lpj_to_grid_func_centre([luh2_dir,'/',luh2_file],1,0));
        luh2_ffrac=luh2_in(:,:,4);
        clear luh2
        gpp=gpp.*luh2_ffrac;
    end
    
    %Get LAI data
    lai=squeeze(lpj_to_grid_func_centre('lai_2001_2014',1,0));
    lai=lai(:,:,1:length(PFTs));
    if ss==1
        lai_1930=squeeze(lpj_to_grid_func_centre('lai_1901_1930',1,0));
        lai_1930=lai_1930(:,:,1:length(PFTs));
    end
    
    %Get disturbance probability data
    distprob_in=squeeze(lpj_to_grid_func_centre('distprob_2001_2014',1,0));
    distprob=distprob_in(:,:,1);
    if ss==1
        temprange=distprob_in(:,:,5);
        WD=distprob_in(:,:,6);
    end
    clear distprob_in
    
    %Get age class data
    ageclass=squeeze(lpj_to_grid_func_centre('ageclass_2014',1,0));
    
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

    %GPP
    if ss==1 || ss==4
        variable='GPP';
        variable_longname='Gross primary productivity';
        units='kg C m-2 yr-1';
        disttype=dtype{ss};
        realisation=dreal{ss};
        modname='LPJ-GUESS';
        year=2014;
        note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';

        write_netcdf_lpjg(gpp,variable,modname,disttype,realisation,year,...
            units,variable_longname,output_dir,note1);
    end
    
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
        units,variable_longname,output_dir,note1,PFTs,'PFT');
    
    if ss==1
        variable='LAI';
        variable_longname='Leaf area index';
        units='m2 (leaf) m-2 (ground)';
        disttype=dtype{ss};
        realisation=dreal{ss};
        modname='LPJ-GUESS';
        year=1930;
        note1='This data is the mean over the period 1961-1990, the year variable in this file is nominal';
        
        write_netcdf_lpjg(lai_1930,variable,modname,disttype,realisation,year,...
            units,variable_longname,output_dir,note1,PFTs,'PFT');
    end
    
    %Disturbance probability
    variable='distprob';
    variable_longname='Disturbance probability';
    units='fractional probability';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';
    
    write_netcdf_lpjg(distprob,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1);
    
    if ss==1
        %Temperature range
        variable='temprange';
        variable_longname='Temperature range';
        units='degrees C';
        disttype=dtype{ss};
        realisation=dreal{ss};
        modname='LPJ-GUESS';
        year=2014;
        note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';

        write_netcdf_lpjg(temprange,variable,modname,disttype,realisation,year,...
            units,variable_longname,output_dir,note1);

        %Wood density
        variable='wooddensity';
        variable_longname='Community mean wood density';
        units='kg DM m-3';
        disttype=dtype{ss};
        realisation=dreal{ss};
        modname='LPJ-GUESS';
        year=2014;
        note1='This data is the mean over the period 2001-2014, the year variable in this file is nominal';

        write_netcdf_lpjg(WD,variable,modname,disttype,realisation,year,...
            units,variable_longname,output_dir,note1);
    end

    %Age class fraction
    variable='age';
    variable_longname='Stand age';
    units='fraction of grid-cell area';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    year=2014;
    note1='This data is for the year 2014';
    
    write_netcdf_lpjg(ageclass,variable,modname,disttype,realisation,year,...
        units,variable_longname,output_dir,note1,ageclasses,'AgeClass');
end
clear ss
