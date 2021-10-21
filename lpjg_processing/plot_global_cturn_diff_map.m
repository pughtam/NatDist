% Plot map of difference in carbon turnover from two different simulations and write processed data out to
% a text file.
%
% Use a tripartite masking system. First mask by any areas where forest is not simulated (based on vegetation biomass). Then mask
% by current forest area. The mask by whether it is in the temperate or boreal biome.
%
% Should be fed LPJ-GUESS tslice files as input.
%
% Dependencies:
% - readmasks_func.m
% - lpj_to_grid_func_centre.m
% - global_grid_area.m
% - wprctile.m (https://uk.mathworks.com/matlabcentral/fileexchange/16920-returns-weighted-percentiles-of-a-sample)
%
% T. Pugh
% 12.01.19

use_cvegmask=false; %Mask by a minimum simulated vegetation biomass density (not recommended for C calculations)
use_fmask=false; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes
farea_opt='luh2'; %Forest area weighting to use, either luh2 or hansen (luh2 recommended)

absplots=false; %Make plots with absolute values
relplots=false; %Make plots with relative values
regstats_gridcell=false; %Calculate regional stats at gridcell level
regstats_sum=true; %Calculate stats based on regional sums

writetxt=false; %Write array to text file
outfile_name='simplemodel_closedcanopy_best_est_100patch_5pClosedCanopyCover_Pugh2019diff.txt';

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';
%lpjg_dir1='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000';
%lpjg_dir2='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000_luh2';

fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='./data/';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file

luh2_dir='./data/';
luh2_file='lu_1700_2015_luh2_aggregate_sum2x2_midpoint_urban_orig_v20_2001_2014';

%--- Read in data ---

% Read cpool and cflux data
cveg1=permute(ncread([lpjg_dir,'/Cveg_LPJ-GUESS_standard_nat_2014.nc'],'Cveg'),[2 1]);
clitter1=permute(ncread([lpjg_dir,'/Clitter_LPJ-GUESS_standard_nat_2014.nc'],'Clitter'),[2 1]);
csoil1=permute(ncread([lpjg_dir,'/Csoil_LPJ-GUESS_standard_nat_2014.nc'],'Csoil'),[2 1]);
npp1=permute(ncread([lpjg_dir,'/NPP_LPJ-GUESS_standard_nat_2014.nc'],'NPP'),[2 1]);
ctau1=cveg1./npp1;


cveg2=permute(ncread([lpjg_dir,'/Cveg_LPJ-GUESS_standard_anthro_2014.nc'],'Cveg'),[2 1]);
clitter2=permute(ncread([lpjg_dir,'/Clitter_LPJ-GUESS_standard_anthro_2014.nc'],'Clitter'),[2 1]);
csoil2=permute(ncread([lpjg_dir,'/Csoil_LPJ-GUESS_standard_anthro_2014.nc'],'Csoil'),[2 1]);
npp2=permute(ncread([lpjg_dir,'/NPP_LPJ-GUESS_standard_anthro_2014.nc'],'NPP'),[2 1]);
ctau2=cveg2./npp2;

[cvegmask,fmask,bmask,ffrac,bmask_temp,bmask_bor]=readmasks_func(use_cvegmask,true,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);

garea=global_grid_area();

farea=ffrac.*garea; %Forest area in m2

%Read LUH2 area fraction data
luh2_in=squeeze(lpj_to_grid_func_centre([luh2_dir,'/',luh2_file],1,0));
luh2_ffrac=luh2_in(:,:,4);
clear luh2

luh2_farea=luh2_ffrac.*garea;

%--- Data processing ---

if strcmp(farea_opt,'luh2')
    %Correct natural simulation by LUH2 prim+sec cover fraction
    cveg1=cveg1.*luh2_ffrac;
    npp1=npp1.*luh2_ffrac;
    csoil1=csoil1.*luh2_ffrac;
    clitter1=clitter1.*luh2_ffrac;
elseif strcmp(farea_opt,'hansen')
    %Convert LUH2 simulation to values per m2 prim+sec cover (for later scaling by Hansen canopy cover fraction)
    cveg2=cveg2./luh2_ffrac;
    npp2=npp2./luh2_ffrac;
    csoil2=csoil2./luh2_ffrac;
    clitter2=clitter2./luh2_ffrac;
else
    error('farea_opt must be set to luh2 or hansen')
end

%Apply masks as specified

if use_cvegmask
    cveg1=cveg1.*cvegmask;
    npp1=npp1.*cvegmask;
    csoil1=csoil1.*cvegmask;
    clitter1=clitter1.*cvegmask;
    ctau1=ctau1.*cvegmask;
    cveg2=cveg2.*cvegmask;
    npp2=npp2.*cvegmask;
    csoil2=csoil2.*cvegmask;
    clitter2=clitter2.*cvegmask;
    ctau2=ctau2.*cvegmask;
end

if use_fmask
    cveg1=cveg1.*fmask;
    npp1=npp1.*fmask;
    csoil1=csoil1.*fmask;
    clitter1=clitter1.*fmask;
    ctau1=ctau1.*fmask;
    cveg2=cveg2.*fmask;
    npp2=npp2.*fmask;
    csoil2=csoil2.*fmask;
    clitter2=clitter2.*fmask;
    ctau2=ctau2.*fmask;
end

if use_bmask
    cveg1=cveg1.*bmask;
    npp1=npp1.*bmask;
    csoil1=csoil1.*bmask;
    clitter1=clitter1.*bmask;
    ctau1=ctau1.*bmask;
    cveg2=cveg2.*bmask;
    npp2=npp2.*bmask;
    csoil2=csoil2.*bmask;
    clitter2=clitter2.*bmask;
    ctau2=ctau2.*bmask;
end

% Absolute difference arrays
diff_cveg=cveg2-cveg1;
diff_ctau=ctau2-ctau1;

% Percentage difference arrays
diffperc_cveg=((cveg2-cveg1)./cveg1)*100;
diffperc_ctau=((ctau2-ctau1)./ctau1)*100;


%--- Make the maps ---
   
if absplots
    % Cveg
    cmin=min(diff_cveg(:));
    cmax=max(diff_cveg(:));
    crange=max(abs([cmin,cmax]));
    cmapmin=-crange-(crange*2)/200;
    
    make_map(ocean_file,cmapmin,crange,diff_cveg);
    title('Cveg absolute difference, Current - Natural')
    
    % Ctau
    cmin=min(diff_ctau(:));
    cmax=max(diff_ctau(:));
    crange=max(abs([cmin,cmax]));
    cmapmin=-crange-(crange*2)/200;
    
    make_map(ocean_file,cmapmin,crange,diff_ctau);
    title('C\tau absolute difference, Current - Natural')
end

if relplots
    % Cveg
    cmin=min(diffperc_cveg(:));
    cmax=max(diffperc_cveg(:));
    crange=min([max(abs([cmin,cmax])),100]);
    cmapmin=-crange-(crange*2)/200;
    
    make_map(ocean_file,cmapmin,crange,diffperc_cveg);
    title('Cveg percentage difference, Current - Natural')
    
    % Ctau
    cmin=min(diffperc_ctau(:));
    cmax=max(diffperc_ctau(:));
    crange=min([max(abs([cmin,cmax])),100]);
    cmapmin=-crange-(crange*2)/200;
    
    make_map(ocean_file,cmapmin,crange,diffperc_ctau);
    title('C\tau percentage difference, Current - Natural')
end


%--- Calculate regional stats ---

if regstats_gridcell
    % Calculate stats across grid cells, weighted by canopy cover area
    if use_fmask==false
        fprintf('Warning: recommended to set use_fmask=true for gridcell-level stats')
    end
    %Initialise weights arrays
    if strcmp(farea_opt,'luh2')
        weights_boreal=luh2_farea(bmask_bor==1);
        weights_temp=luh2_farea(bmask_temp==1);
        weights_boreal_litt=luh2_farea(bmask_bor==1);
        weights_temp_litt=luh2_farea(bmask_temp==1);
    elseif strcmp(farea_opt,'hansen')
        weights_boreal=farea(bmask_bor==1);
        weights_temp=farea(bmask_temp==1);
        weights_boreal_litt=farea(bmask_bor==1);
        weights_temp_litt=farea(bmask_temp==1);
    end
    
    diffperc_ctau_boreal=diffperc_ctau(bmask_bor==1);
    weights_boreal(isnan(diffperc_ctau_boreal))=[];
    diffperc_ctau_boreal(isnan(diffperc_ctau_boreal))=[];
    
    diffperc_ctau_temp=diffperc_ctau(bmask_temp==1);
    weights_temp(isnan(diffperc_ctau_temp))=[];
    diffperc_ctau_temp(isnan(diffperc_ctau_temp))=[];
    
    diffperc_ctau_median_boreal=wmedian(diffperc_ctau_boreal,weights_boreal);
    diffperc_ctau_median_temp=wmedian(diffperc_ctau_temp,weights_temp);
    
    diffperc_ctau_mean_boreal=wmean(diffperc_ctau_boreal,weights_boreal);
    diffperc_ctau_mean_temp=wmean(diffperc_ctau_temp,weights_temp);
    
    diffperc_ctau_10perc_boreal=wprctile(diffperc_ctau_boreal,10,weights_boreal);
    diffperc_ctau_90perc_boreal=wprctile(diffperc_ctau_boreal,90,weights_boreal);
    diffperc_ctau_10perc_temp=wprctile(diffperc_ctau_temp,10,weights_temp);
    diffperc_ctau_90perc_temp=wprctile(diffperc_ctau_temp,90,weights_temp);
    
    
    diffperc_clitter_boreal=((clitter2(bmask_bor==1)-clitter1(bmask_bor==1))./clitter1(bmask_bor==1))*100;
    weights_boreal_litt(isnan(diffperc_clitter_boreal))=[];
    diffperc_clitter_boreal(isnan(diffperc_clitter_boreal))=[];
    
    diffperc_clitter_temp=((clitter2(bmask_temp==1)-clitter1(bmask_temp==1))./clitter1(bmask_temp==1))*100;
    weights_temp_litt(isnan(diffperc_clitter_temp))=[];
    diffperc_clitter_temp(isnan(diffperc_clitter_temp))=[];
    
    diffperc_clitter_median_boreal=wmedian(diffperc_clitter_boreal,weights_boreal_litt);
    diffperc_clitter_median_temp=wmedian(diffperc_clitter_temp,weights_temp_litt);
    
    diffperc_clitter_mean_boreal=wmean(diffperc_clitter_boreal,weights_boreal_litt);
    diffperc_clitter_mean_temp=wmean(diffperc_clitter_temp,weights_temp_litt);
    
    diffperc_clitter_10perc_boreal=wprctile(diffperc_clitter_boreal,10,weights_boreal_litt);
    diffperc_clitter_90perc_boreal=wprctile(diffperc_clitter_boreal,90,weights_boreal_litt);
    diffperc_clitter_10perc_temp=wprctile(diffperc_clitter_temp,10,weights_temp_litt);
    diffperc_clitter_90perc_temp=wprctile(diffperc_clitter_temp,90,weights_temp_litt);
end

if regstats_sum
    % Calculate mean biome-level stats
    if use_fmask==true
        fprintf('Warning: recommended to set use_fmask=false for biome-level stats')
    end
    
    if strcmp(farea_opt,'luh2')
        %Values are already per gridcell, just need multiplying by gridcell area
        cveg1_area=cveg1.*garea;
        cveg2_area=cveg2.*garea;
        csoil1_area=csoil1.*garea;
        csoil2_area=csoil2.*garea;
        clitter1_area=clitter1.*garea;
        clitter2_area=clitter2.*garea;
        npp1_area=npp1.*garea;
        npp2_area=npp2.*garea;
    elseif strcmp(farea_opt,'hansen')
        %Values are per prim+sec (natural) area and need correcting by forest cover fraction (here canopy cover fraction)
        cveg1_area=cveg1.*farea;
        cveg2_area=cveg2.*farea;
        csoil1_area=csoil1.*farea;
        csoil2_area=csoil2.*farea;
        clitter1_area=clitter1.*farea;
        clitter2_area=clitter2.*farea;
        npp1_area=npp1.*farea;
        npp2_area=npp2.*farea;
    end
    
    cveg1_boreal=cveg1_area(bmask_bor==1);
    cveg2_boreal=cveg2_area(bmask_bor==1);
    npp1_boreal=npp1_area(bmask_bor==1);
    npp2_boreal=npp2_area(bmask_bor==1);
    csoil1_boreal=csoil1_area(bmask_bor==1);
    csoil2_boreal=csoil2_area(bmask_bor==1);
    clitter1_boreal=clitter1_area(bmask_bor==1);
    clitter2_boreal=clitter2_area(bmask_bor==1);
    
    ctau1_boreal=nansum(cveg1_boreal)/nansum(npp1_boreal);
    ctau2_boreal=nansum(cveg2_boreal)/nansum(npp2_boreal);
    ctau_diff_allboreal_perc=((ctau2_boreal-ctau1_boreal)/ctau1_boreal)*100;
    ctaueco1_boreal=(nansum(csoil1_boreal)+nansum(clitter1_boreal)+nansum(cveg1_boreal))/nansum(npp1_boreal);
    ctaueco2_boreal=(nansum(csoil2_boreal)+nansum(clitter2_boreal)+nansum(cveg2_boreal))/nansum(npp2_boreal);
    ctaueco_diff_allboreal_perc=((ctaueco2_boreal-ctaueco1_boreal)/ctaueco1_boreal)*100;
    cveg1_sum_boreal=nansum(cveg1_boreal)/1e12; %In Pg C
    cveg2_sum_boreal=nansum(cveg2_boreal)/1e12; %In Pg C
    csoil1_sum_boreal=nansum(csoil1_boreal)/1e12; %In Pg C
    csoil2_sum_boreal=nansum(csoil2_boreal)/1e12; %In Pg C
    clitter1_sum_boreal=nansum(clitter1_boreal)/1e12; %In Pg C
    clitter2_sum_boreal=nansum(clitter2_boreal)/1e12; %In Pg C
    
    cveg1_temp=cveg1_area(bmask_temp==1);
    cveg2_temp=cveg2_area(bmask_temp==1);
    npp1_temp=npp1_area(bmask_temp==1);
    npp2_temp=npp2_area(bmask_temp==1);
    csoil1_temp=csoil1_area(bmask_temp==1);
    csoil2_temp=csoil2_area(bmask_temp==1);
    clitter1_temp=clitter1_area(bmask_temp==1);
    clitter2_temp=clitter2_area(bmask_temp==1);
    weights_temp=ffrac(bmask_temp==1);
    
    ctau1_temp=nansum(cveg1_temp)/nansum(npp1_temp);
    ctau2_temp=nansum(cveg2_temp)/nansum(npp2_temp);
    ctau_diff_alltemp_perc=((ctau2_temp-ctau1_temp)/ctau1_temp)*100;
    ctaueco1_temp=(nansum(csoil1_temp)+nansum(clitter1_temp)+nansum(cveg1_temp))/nansum(npp1_temp);
    ctaueco2_temp=(nansum(csoil2_temp)+nansum(clitter2_temp)+nansum(cveg2_temp))/nansum(npp2_temp);
    ctaueco_diff_alltemp_perc=((ctaueco2_temp-ctaueco1_temp)/ctaueco1_temp)*100;
    cveg1_sum_temp=nansum(cveg1_temp)/1e12; %In Pg C
    cveg2_sum_temp=nansum(cveg2_temp)/1e12; %In Pg C
    csoil1_sum_temp=nansum(csoil1_temp)/1e12; %In Pg C
    csoil2_sum_temp=nansum(csoil2_temp)/1e12; %In Pg C
    clitter1_sum_temp=nansum(clitter1_temp)/1e12; %In Pg C
    clitter2_sum_temp=nansum(clitter2_temp)/1e12; %In Pg C
end


%--- Write out to text file ---

if writetxt
    plotarray_out=plotarray;
    plotarray_out(isnan(plotarray))=-9999;
    plotarray_out=int32(plotarray_out);
    
    fid=fopen(outfile_name,'w');
    for yy=360:-1:1
        fprintf(fid,repmat('%7d',1,720),plotarray_out(yy,:));
        fprintf(fid,'\n');
    end
    clear yy
    fclose(fid);
end


%--- Plotting function ---
function make_map(ocean_file,cmapmin,crange,plotarray)

%Read mask for ocean areas
    load(ocean_file,'esa_05');
    oceanm=NaN(720,360);
    oceanm(esa_05'>200 & esa_05'<220)=cmapmin;
    oceanm=oceanm';
    
    lons=-180:0.5:179.5;
    lats=-90:0.5:89.5;
    
    figure
    %cmap=colormap(redblue(200));
    cmap=cbrewer('div','RdYlBu',100);
    cmap=[0.9 0.9 0.9; cmap];
    
    colormap(cmap)
    axesm('MapProjection','robinson','MapLatLimit',[23 80])
    hold on
    l1=pcolorm(lats,lons,oceanm);
    set(l1,'linestyle','none')
    p1=pcolorm(lats,lons,plotarray);
    set(p1,'linestyle','none')
    axis tight
    caxis([cmapmin crange])
    c1=colorbar;
    set(c1,'FontSize',12,'FontWeight','Bold')
    set(c1,'Limits',[-crange crange])
end
