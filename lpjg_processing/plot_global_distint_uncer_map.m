% Plot map of uncertainty in disturbance rotation periods, showing the absolute range of uncertainty divided by the best
% estimate.
%
% Use a tripartite masking system. First mask by any areas where forest is not simulated (based on vegetation biomass). Then mask
% by current forest area. The mask by whether it is in the temperate or boreal biome.
%
% Dependencies:
% - readmasks_func.m
%
% T. Pugh
% 23.01.20

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

limitscale=false; %Cap disturbance interval at 1000 years
dimplot=1; %Column number in input file containing the disturbance rate

makeplot=true; %Make a plot
writetxt=false; %Write array to text file
outfile_name='simplemodel_best_est_10pCanopyCover_uncerfrac.txt';

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';
lpjg_file_base='distprob_LPJ-GUESS_standard_nat_2014.nc';
lpjg_file_low='distprob_LPJ-GUESS_low_nat_2014.nc';
lpjg_file_high='distprob_LPJ-GUESS_high_nat_2014.nc';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='./data/';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

% Read disturbance interval
dist_base=permute(ncread([lpjg_dir,'/',lpjg_file_base],'distprob'),[2 1]); 
dist_base(dist_base==0)=NaN;

dist_low=permute(ncread([lpjg_dir,'/',lpjg_file_low],'distprob'),[2 1]); 
dist_low(dist_low==0)=NaN;

dist_high=permute(ncread([lpjg_dir,'/',lpjg_file_high],'distprob'),[2 1]); 
dist_high(dist_high==0)=NaN;

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


%--- Data processing ---

%Create an array for plotting the disturbance interval

distint_base=1./dist_base(:,:,dimplot);
distint_low=1./dist_low(:,:,dimplot);
distint_high=1./dist_high(:,:,dimplot);

if limitscale
    distint_base(distint_base>1000)=1000;
    distint_low(distint_low>1000)=1000;
    distint_high(distint_high>1000)=1000;
end

plotarray=abs((distint_low-distint_high)./distint_base);

if use_cvegmask
    plotarray=plotarray.*cvegmask;
end

if use_fmask
    plotarray=plotarray.*fmask;
end

if use_bmask
    plotarray=plotarray.*bmask;
end

cmax=0.3;
cmin=0;

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
    cmap=colormap(parula(200));
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
    
end


%--- Write out to text file ---

if writetxt
    plotarray_out=plotarray;
    plotarray_out(isnan(plotarray))=-9999;
    
    fid=fopen(outfile_name,'w');
    for yy=360:-1:1
        fprintf(fid,repmat('%7.3f',1,720),plotarray_out(yy,:));
        fprintf(fid,'\n');
    end
    clear yy
    fclose(fid);
end
