% Plot map of difference between disturbance rotation periods from two different simulations and write processed data out to
% a text file.
%
% Use a tripartite masking system. First mask by any areas where forest is not simulated (based on vegetation biomass). Then mask
% by current forest area. The mask by whether it is in the temperate or boreal biome.
%
% Should be fed LPJ-GUESS tslice files as input.
%
% Dependencies:
% - readmasks_func.m
%
% T. Pugh
% 12.01.19

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=true; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

logscale=false; %Plot map using quasi-log colour scale
limitscale=false; %Cap disturbance interval at 1000 years
dimplot=1; %Column number in input file containing the disturbance rate

makeplot=true; %Make a plot
writetxt=true; %Write array to text file
outfile_name='simplemodel_best_est_5pClosedCanopyCover_diffClosedCanopy.txt';

lpjg_dir1='/Users/pughtam/LPJG/disturbance_prognostic_runs/simplemodel_best_est';
lpjg_dir2='/Users/pughtam/LPJG/disturbance_prognostic_runs/simplemodel_closedcanopy_best_est';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

% Read disturbance interval
dist1=squeeze(lpj_to_grid_func_centre([lpjg_dir1,'/distprob_2001_2014'],1,0));
dist1(dist1==0)=NaN;

dist2=squeeze(lpj_to_grid_func_centre([lpjg_dir2,'/distprob_2001_2014'],1,0));
dist2(dist2==0)=NaN;

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir1,fmask_dir,fmask_file,bmask_dir);


%--- Data processing ---

%Create an array for plotting the disturbance interval

%Calculate the difference between the two disturbance interval sets
distint1=1./dist1(:,:,dimplot);
distint2=1./dist2(:,:,dimplot);

if limitscale
    distint1(distint1>1000)=1000;
    distint2(distint2>1000)=1000;
end

plotarray=distint2-distint1;

if use_cvegmask
    plotarray=plotarray.*cvegmask;
end

if use_fmask
    plotarray=plotarray.*fmask;
end

if use_bmask
    plotarray=plotarray.*bmask;
end

%Create a quasi-logarithmic scale that displays both positive and negative changes
if logscale
    plotarray_log=log10(abs(plotarray));
    plotarray_log(plotarray_log<0)=0;
    plotarray_log(plotarray<0)=-plotarray_log(plotarray<0);
    plotarray=plotarray_log;
    clear plotarray_log
    
    cmax=3;
    cmin=-3;
else
    cmax=1000;
    cmin=-1000;
end

crange=cmax-cmin;
cmapmin=cmin-crange/200;


%--- Make the map ---

if makeplot
    
    %Read mask for ocean areas
    load(ocean_file);
    oceanm=NaN(720,360);
    oceanm(esa_05'>200 & esa_05'<220)=cmapmin;
    oceanm=oceanm';
    
    lons=-180:0.5:179.5;
    lats=-90:0.5:89.5;
    
    figure
    cmap=colormap(redblue(200));
    cmap=[0.9 0.9 0.9; cmap];
    
    colormap(cmap)
    axesm('MapProjection','robinson','MapLatLimit',[23 80])
    hold on
    l1=pcolorm(lats,lons,oceanm);
    set(l1,'linestyle','none')
    p1=pcolorm(lats,lons,plotarray);
    set(p1,'linestyle','none')
    axis tight
    caxis([cmapmin cmax])
    c1=colorbar;
    set(c1,'FontSize',12,'FontWeight','Bold')
    set(c1,'Limits',[cmin cmax])
    if logscale
        set(c1,'Ticks',cmin:1:cmax);
        set(c1,'TickLabels',[-1000 -100 -10 0 10 100 1000])
    end
    
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