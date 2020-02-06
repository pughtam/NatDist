% Plot map of disturbance rotation periods and write processed data out to a text file
%
% Use a tripartite masking system. First mask by any areas where forest is not simulated (based on vegetation biomass). Then mask
% by current forest area. The mask by whether it is in the temperate or boreal biome.
%
% Should be fed an LPJ-GUESS tslice file as input.
%
% Dependencies:
% - readmasks_func.m
%
% T. Pugh
% 12.01.19

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

distvar=true; %Plot as disturbance return interval (true) or raw variable (false)
logscale=true; %Plot map using log colour scale
limitscale=false; %Cap colour scale at 1000 years
dimplot=1; %Column number in input file containing the disturbance rate

makeplot=true; %Make a plot
writetxt=false; %Write array to text file
writenetcdf=true; %Write array to netcdf file
outfile_name='simplemodel_best_est_100patch_10pCanopyCover_nolimit';

lpjg_dir='/Users/pughtam/LPJG/disturbance_prognostic_runs/simplemodel_best_est_100patch';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

% Read disturbance interval
dist=squeeze(lpj_to_grid_func_centre([lpjg_dir,'/distprob_2001_2014'],1,0));
dist(dist==0)=NaN;

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


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
    
    fid=fopen([outfile_name,'.txt'],'w');
    for yy=360:-1:1
        fprintf(fid,repmat('%7d',1,720),plotarray_out(yy,:));
        fprintf(fid,'\n');
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
    
    minlat=-89.75;
    maxlat=89.75;
    minlon=-179.75;
    maxlon=179.75;
    gridspace=0.5;
    
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