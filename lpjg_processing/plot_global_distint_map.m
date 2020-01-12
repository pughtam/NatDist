% Plot map of disturbance rotation periods and write processed data out to a text file
%
% Use a dual masking system. First mask by any areas where forest is not simulated (based on vegetation biomass). Then mask
% by current forest area.
%
% Should be fed an LPJ-GUESS tslice file as input.
%
% T. Pugh
% 12.01.19

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=true; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

distvar=true; %Plot as disturbance return interval (true) or raw variable (false)
logscale=true; %Plot map using log colour scale
limitscale=false; %Cap colour scale at 1000 years
dimplot=1; %Column number in input file containing the disturbance rate

makeplot=false; %Make a plot
writetxt=true; %Write array to text file
outfile_name='simplemodel_closedcanopy_low_est_5pClosedCanopyCover_nolimit.txt';

lpjg_dir='/Users/pughtam/LPJG/disturbance_prognostic_runs/simplemodel_closedcanopy_low_est';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

% Read disturbance interval
dist=squeeze(lpj_to_grid_func_centre([lpjg_dir,'/distprob_2001_2014'],1,0));
dist(dist==0)=NaN;

% Create mask to exclude grid cells where the vegetation C mass does not meet a threshold of 1 kg C m-2
if use_cvegmask
    cpool=lpj_to_grid_func_centre([lpjg_dir,'/cpool_2001_2014'],1,0);
    cveg=squeeze(cpool(:,:,1));
    cvegmask=NaN(size(cveg));
    cvegmask(cveg>1)=1;
end

% Create mask to exclude grid cells where at least 10% canopy is not reached
if use_fmask
    if ccmask
        fmask_var='forested_50_percent'; %Name of variable in fmask file to use
        ffrac=ncread([fmask_dir,'/',fmask_file],fmask_var)';
        fmask=NaN(size(ffrac));
        fmask(ffrac>=5)=1;
    else
        fmask_var='canopy_cover_frac'; %Name of variable in fmask file to use
        ffrac=ncread([fmask_dir,'/',fmask_file],fmask_var)';
        fmask=NaN(size(ffrac));
        fmask(ffrac>=10)=1;
    end
end

%Read in the biome mask and format to 0.5 x 0.5 degrees
if use_bmask
    bmask_temp=flipud(geotiffread([bmask_dir,'/temperate_biome_025degree.tif']));
    bmask_bor=flipud(geotiffread([bmask_dir,'/boreal_biome_025degree.tif']));
    bmask_025=zeros(size(bmask_temp));
    bmask_025(bmask_temp==1)=1;
    bmask_025(bmask_bor==1)=1;
    clear bmask_bor bmask_temp
    bmask=NaN(360,720);
    for xx=1:720
        for yy=1:360
            xx_s=(xx*2)-1;
            xx_e=xx*2;
            yy_s=(yy*2)-1;
            yy_e=yy*2;
            temp=bmask_025(yy_s:yy_e,xx_s:xx_e);
            bmask(yy,xx)=mode(temp(:));
        end
    end
    bmask(bmask==0)=NaN;
    clear xx yy xx_s xx_e yy_s yy_e    
    clear bmask_025
end


%--- Data processing ---

%Create an array for plotting the disturbance interval
if logscale && distvar
    plotarray=log10(1./squeeze(dist(:,:,dimplot)));
else
    if distvar
        plotarray=1./squeeze(dist(:,:,dimplot));
    else
        plotarray=squeeze(dist(:,:,dimplot));
    end
end

if use_cvegmask
    plotarray=plotarray.*cvegmask;
end

if use_fmask
    plotarray=plotarray.*fmask;
end

if use_bmask
    plotarray=plotarray.*bmask;
end

% Set caxis ranges
if logscale
    cmin=log10(10);
    if limitscale
        cmax=log10(1000);
    else
        cmax=log10(4000);
    end
    crange=cmax-cmin;
    cmapmin=cmin-log10(1);
else
    if distvar
        cmin=10;
        if limitscale
            cmax=4000;
        else
            cmax=1000;
        end
    else
        cmin=0;
        cmax=1;
    end
    crange=cmax-cmin;
    cmapmin=-crange/200;
end


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
    if logscale
        set(c1,'Ticks',cmin:0.5:cmax);
        set(c1,'TickLabels',round(10.^(get(c1,'Ticks'))))
    end
    
end


%--- Write out to text file ---

if writetxt
    
    if logscale
        plotarray_out=10.^(plotarray);
    else
        plotarray_out=plotarray;
    end
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