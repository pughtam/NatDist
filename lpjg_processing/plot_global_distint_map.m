% Plot map of disturbance rotation periods and write processed data out to a text file
%
% Use a tripartite masking system. First mask by any areas where forest is not simulated (based on vegetation biomass). Then mask
% by current forest area. The mask by whether it is in the temperate or boreal biome.
%
% Should be fed an LPJ-GUESS tslice file as input.
%
% Dependencies:
% - readmasks_func.m
% - lpj_to_grid_func_centre.m
%
% T. Pugh
% 12.01.19

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=true; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

distvar=true; %Plot as disturbance return interval (true) or raw variable (false)
logscale=false; %Plot map using log colour scale
limitscale=true; %Cap colour scale at 1000 years
dimplot=1; %Column number in input file containing the disturbance rate

makeplot=true; %Make a plot
readnetcdf=false; %Read from a netcdf file, otherwise from an LPJ-GUESS output file
writetxt=false; %Write array to text file
writenetcdf=false; %Write array to netcdf file
output1deg=false; %Write netcdf at 1 x 1 degree aggregation, instead of 0.5 x 0.5
%outfile_name='best_est_adjparam_latosa4000_closedcan_20patch_5pClosedCanopyCover_1deg';
outfile_name='best_est_adjparam_latosa4000_20patch_10pCanopyCover';
makeregionstats=true; %Make stats at regional level

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';
%lpjg_file='distprob_LPJ-GUESS_standard_nat_2014.nc';
lpjg_file='distprob_LPJ-GUESS_standard_natcc_2014.nc';
netcdf_file='/Users/pughtam/Documents/GAP_and_other_work/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc';
netcdf_varname='tauO';

fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='./data/';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

% Read disturbance interval
if readnetcdf
    % Read disturbance interval from netcdf
    distint_1deg=ncread(netcdf_file,netcdf_varname)';
    % Resample to 0.5 degree
    distint=NaN(360,720);
    for xx=1:720
        for yy=1:360
            xx_s=ceil(xx/2);
            yy_s=ceil(yy/2);
            distint(yy,xx)=distint_1deg(yy_s,xx_s);
        end
    end
    clear xx yy xx_s yy_s distint_1deg
    dist=1./distint;
    clear distint
else
    dist=permute(ncread([lpjg_dir,'/',lpjg_file],'distprob'),[2 1]); 
    dist(dist==0)=NaN;
end

[cvegmask,fmask,bmask,ffrac]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


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
            cmax=1000;
        else
            cmax=4000;
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

%--- Aggregate to one degree for output? ---

if output1deg
    plotarray_1deg=NaN(180,360);
    for xx=1:360
        for yy=1:180
            xx_s=(xx*2)-1;
            xx_e=xx*2;
            yy_s=(yy*2)-1;
            yy_e=yy*2;
            temp=plotarray(yy_s:yy_e,xx_s:xx_e);
            plotarray_1deg(yy,xx)=nanmean(temp(:));
        end
    end
    clear xx yy xx_s xx_e yy_s yy_e temp
    plotarray=plotarray_1deg;
    clear plotarray_1deg
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
    
    fid=fopen([outfile_name,'.txt'],'w');
    if output1deg
        for yy=180:-1:1
            fprintf(fid,repmat('%7d',1,360),plotarray_out(yy,:));
            fprintf(fid,'\n');
        end
    else
        for yy=360:-1:1
            fprintf(fid,repmat('%7d',1,720),plotarray_out(yy,:));
            fprintf(fid,'\n');
        end
    end
    clear yy
    fclose(fid);
end

%--- Write out to netcdf file ---

if writenetcdf
    
    fprintf('Creating netcdf file\n')
    
    if logscale
        plotarray_out=10.^(plotarray);
    else
        plotarray_out=plotarray;
    end
    
    if output1deg
        minlat=-89.5;
        maxlat=89.5;
        minlon=-179.5;
        maxlon=179.5;
        gridspace=1.0;
    else
        minlat=-89.75;
        maxlat=89.75;
        minlon=-179.75;
        maxlon=179.75;
        gridspace=0.5;
    end
    
    latgrid=minlat:gridspace:maxlat;
    longrid=minlon:gridspace:maxlon;
    nlat=length(latgrid);
    nlon=length(longrid);
    
    outfile=[outfile_name,'.nc'];
    
    nccreate(outfile,'latitude','Dimensions',{'latitude',nlat})
    ncwrite(outfile,'latitude',latgrid)
    ncwriteatt(outfile,'latitude','units','degrees_north');
    nccreate(outfile,'longitude','Dimensions',{'longitude',nlon})
    ncwrite(outfile,'longitude',longrid)
    ncwriteatt(outfile,'longitude','units','degrees_east');
    
    nccreate(outfile,'tau','Dimensions',{'longitude','latitude'},'DeflateLevel',9)
    
    ncwrite(outfile,'tau',plotarray_out')
    ncwriteatt(outfile,'tau','longname','Modelled disturbance rotation time')
    ncwriteatt(outfile,'tau','units','years')
    
    ncwriteatt(outfile,'/','Institution','University of Birmingham, UK');
    ncwriteatt(outfile,'/','Contact','Thomas Pugh, t.a.m.pugh@bham.ac.uk');
    ncwriteatt(outfile,'/','Version',['Version 1: ',date]);
    
end


%--- Make some regional statistics ---

if makeregionstats
    % Read and prepare the region mask
    bmask_temp=flipud(geotiffread([bmask_dir,'/temperate_biome_025degree.tif']));
    bmask_bor=flipud(geotiffread([bmask_dir,'/boreal_biome_025degree.tif']));
    bmask_025=zeros(size(bmask_temp));
    bmask_025(bmask_temp==1)=1;
    bmask_025(bmask_bor==1)=2;
    rmask=NaN(180,360);
    for xx=1:720
        for yy=1:360
            xx_s=(xx*2)-1;
            xx_e=xx*2;
            yy_s=(yy*2)-1;
            yy_e=yy*2;
            temp=bmask_025(yy_s:yy_e,xx_s:xx_e);
            rmask(yy,xx)=mode(temp(:));
        end
    end
    rmask(rmask==0)=NaN;
    clear xx yy xx_s xx_e yy_s yy_e
    clear bmask_025
    clear temp
    nregion=2;
    regions={'Boreal forest','Temperate forest'};
    regions_short={'boreal','temperate'};
    
    % Get forest area (for weighting)
    garea=global_grid_area();
    farea=ffrac.*garea;
    
    % Calculate stats
    plotarray_median_temp=wmedian(plotarray(rmask==1),farea(rmask==1));
    plotarray_median_boreal=wmedian(plotarray(rmask==2),farea(rmask==2));
    
end
