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

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

absplots=false; %Make plots with absolute values
relplots=true; %Make plots with relative values
regstats=true; %Calculate regional stats

writetxt=false; %Write array to text file
outfile_name='simplemodel_closedcanopy_best_est_100patch_5pClosedCanopyCover_Pugh2019diff.txt';

lpjg_dir1='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000';
lpjg_dir2='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000_luh2';

fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

% Read cpool and cflux data
cpool1=squeeze(lpj_to_grid_func_centre([lpjg_dir1,'/cpool_2001_2014'],1,0));
cveg1=cpool1(:,:,1);
csoil1=sum(cpool1(:,:,2:3),3);
clear cpool1

cflux1=squeeze(lpj_to_grid_func_centre([lpjg_dir1,'/cflux_2001_2014'],1,0));
crepro1=cflux1(:,:,2);
npp1=-cflux1(:,:,1);
ctau1=cveg1./npp1;
clear cflux1


cpool2=squeeze(lpj_to_grid_func_centre([lpjg_dir2,'/cpool_2001_2014'],1,0));
cveg2=cpool2(:,:,1);
csoil2=sum(cpool2(:,:,2:3),3);
clear cpool2
    
cflux2=squeeze(lpj_to_grid_func_centre([lpjg_dir2,'/cflux_2001_2014'],1,0));
crepro2=cflux2(:,:,2);
npp2=-cflux2(:,:,1);
ctau2=cveg2./npp2;
clear cflux2

[cvegmask,fmask,bmask,ffrac,bmask_temp,bmask_bor]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir1,fmask_dir,fmask_file,bmask_dir);

garea=global_grid_area();

farea=ffrac.*garea; %Forest area in m2

%--- Data processing ---

%Apply masks as specified

if use_cvegmask
    cveg1=cveg1.*cvegmask;
    npp1=npp1.*cvegmask;
    csoil1=csoil1.*cvegmask;
    ctau1=ctau1.*cvegmask;
    cveg2=cveg2.*cvegmask;
    npp2=npp2.*cvegmask;
    csoil2=csoil2.*cvegmask;
    ctau2=ctau2.*cvegmask;
end

if use_fmask
    cveg1=cveg1.*fmask;
    npp1=npp1.*fmask;
    csoil1=csoil1.*fmask;
    ctau1=ctau1.*fmask;
    cveg2=cveg2.*fmask;
    npp2=npp2.*fmask;
    csoil2=csoil2.*fmask;
    ctau2=ctau2.*fmask;
end

if use_bmask
    cveg1=cveg1.*bmask;
    npp1=npp1.*bmask;
    csoil1=csoil1.*bmask;
    ctau1=ctau1.*bmask;
    cveg2=cveg2.*bmask;
    npp2=npp2.*bmask;
    csoil2=csoil2.*bmask;
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

% Calculate stats across grid cells, weighted by forest cover area
diffperc_ctau_boreal=diffperc_ctau(bmask_bor==1);
weights_boreal=ffrac(bmask_bor==1);
weights_boreal(isnan(diffperc_ctau_boreal))=[];
diffperc_ctau_boreal(isnan(diffperc_ctau_boreal))=[];

diffperc_ctau_temp=diffperc_ctau(bmask_temp==1);
weights_temp=ffrac(bmask_temp==1);
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

% Calculate mean biome-level stats
cveg1_area=cveg1.*farea;
cveg2_area=cveg2.*farea;
csoil1_area=csoil1.*farea;
csoil2_area=csoil2.*farea;
npp1_area=npp1.*farea;
npp2_area=npp2.*farea;

cveg1_boreal=cveg1_area(bmask_bor==1);
cveg2_boreal=cveg2_area(bmask_bor==1);
npp1_boreal=npp1_area(bmask_bor==1);
npp2_boreal=npp2_area(bmask_bor==1);
csoil1_boreal=csoil1_area(bmask_bor==1);
csoil2_boreal=csoil2_area(bmask_bor==1);

ctau1_boreal=nansum(cveg1_boreal)/nansum(npp1_boreal);
ctau2_boreal=nansum(cveg2_boreal)/nansum(npp2_boreal);
ctau_diff_allboreal_perc=((ctau2_boreal-ctau1_boreal)/ctau1_boreal)*100;
ctaueco1_boreal=nansum(csoil1_boreal)/nansum(npp1_boreal);
ctaueco2_boreal=nansum(csoil2_boreal)/nansum(npp2_boreal);
ctaueco_diff_allboreal_perc=((ctaueco2_boreal-ctaueco1_boreal)/ctaueco1_boreal)*100;
cveg1_sum_boreal=nansum(cveg1_boreal);
cveg2_sum_boreal=nansum(cveg2_boreal);

cveg1_temp=cveg1_area(bmask_temp==1);
cveg2_temp=cveg2_area(bmask_temp==1);
npp1_temp=npp1_area(bmask_temp==1);
npp2_temp=npp2_area(bmask_temp==1);
csoil1_temp=csoil1_area(bmask_temp==1);
csoil2_temp=csoil2_area(bmask_temp==1);
weights_temp=ffrac(bmask_temp==1);

ctau1_temp=nansum(cveg1_temp)/nansum(npp1_temp);
ctau2_temp=nansum(cveg2_temp)/nansum(npp2_temp);
ctau_diff_alltemp_perc=((ctau2_temp-ctau1_temp)/ctau1_temp)*100;
ctaueco1_temp=nansum(csoil1_temp)/nansum(npp1_temp);
ctaueco2_temp=nansum(csoil2_temp)/nansum(npp2_temp);
ctaueco_diff_alltemp_perc=((ctaueco2_temp-ctaueco1_temp)/ctaueco1_temp)*100;
cveg1_sum_temp=nansum(cveg1_temp);
cveg2_sum_temp=nansum(cveg2_temp);


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
