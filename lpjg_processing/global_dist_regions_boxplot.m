% Plot boxplot of disturbance rotation periods by region.
% Regions are based on an aggregation of the GFAD regions.
%
% Dependencies:
% - readmasks_func.m
% - gfad_regions.m
%
%T. Pugh
%24.01.20

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=true; %Use a closed-canopy forest mask (if use_fmask=true and readlpjgdata=true)
use_bmask=true; %Mask by temperate/boreal biomes

limitscale=false; %Cap disturbance interval at 1000 years
dimplot=1;

makeregionmap=false; %Make a map of where the regions are located.

readlpjgdata=false; %Read data from LPJ-GUESS (true) or form a netcdf file (false)

%If reading output from LPJ-GUESS must specify folder (readlpjgdata=true)
lpjg_dir='/Users/pughtam/LPJG/disturbance_prognostic_runs/simplemodel_closedcanopy_best_est_100patch';
%If reading output from netcdf file must specify file location and variable name (readlpjgdata=false)
netcdf_file='/Users/pughtam/Documents/GAP_and_other_work/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc';
netcdf_varname='tauO';
%Ancillary data files
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

if readlpjgdata
    % Read disturbance interval from LPJ-GUESS output
    dist=lpj_to_grid_func_centre([lpjg_dir,'/distprob_2001_2014'],1,0);
    dist(dist==0)=NaN;
    distint=1./dist(:,:,dimplot);
else
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
end

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


gfad_filepath_stan='/Users/pughtam/data/GFAD_V1-1/GFAD_V1-1.nc';
onedeg=false;
[new_region,regionnames]=gfad_regions(gfad_filepath_stan,onedeg);
new_region=flipud(new_region');

%--- Data processing ---

%Create an array for plotting the disturbance interval

if limitscale
    distint(distint>1000)=1000;
end

if use_cvegmask && readlpjgdata
    distint=distint.*cvegmask;
    new_region=new_region.*cvegmask;
end

if use_fmask
    distint=distint.*fmask;
    new_region=new_region.*fmask;
end

if use_bmask
    distint=distint.*bmask;
    new_region=new_region.*bmask;
end

%Modify regions to merge some together
new_region(new_region==19)=1; %Remove Kazakstan
new_region(new_region==18)=1; %Remove Mongolia
new_region(new_region==17)=1; %Remove New Zealand
new_region(new_region==11)=10; %Merge Mid and South China
new_region(new_region>10)=new_region(new_region>10)-1;
new_region(new_region==0)=NaN;
regionnames(17:19)=[];
regionnames(11)=[];

cmin=1;
cmax=max(new_region(:));

crange=cmax-cmin;
cmapmin=cmin-crange/200;


%--- Make a map of the region locations ---

if makeregionmap
    %Read mask for ocean areas
    load(ocean_file);
    oceanm=NaN(720,360);
    oceanm(esa_05'>200 & esa_05'<220)=cmapmin;
    oceanm=oceanm';
    
    lons=-180:0.5:179.5;
    lats=-90:0.5:89.5;
    
    figure
    cmap=cbrewer('qual','Paired',cmax+1);
    cmap=[0.9 0.9 0.9; cmap];
    
    colormap(cmap)
    axesm('MapProjection','robinson','MapLatLimit',[23 80])
    hold on
    l1=pcolorm(lats,lons,oceanm);
    set(l1,'linestyle','none')
    p1=pcolorm(lats,lons,new_region);
    set(p1,'linestyle','none')
    axis tight
    caxis([-1 cmax+0.5])
    c1=colorbar;
    set(c1,'FontSize',12,'FontWeight','Bold')
    set(c1,'Limits',[cmin cmax+0.5])
    set(c1,'Ticks',cmin:cmax)
    set(c1,'TickLabels',regionnames)
end


%--- Make a boxplot (excluding the "rest of the world" category) ---

figure
boxplot(distint(new_region~=1),new_region(new_region~=1),'notch','on','plotstyle','compact')
ylabel('\tau (years)')
set(gca,'XTick',1:cmax-1,'XTickLabel',regionnames(2:cmax))
set(gca,'XTickLabelRotation',300)
set(gca,'YScale','log')
set(gca,'YLim',[100 4000])
set(gca,'YTickLabel',{'100','1000'})
%t4=text(0.7,27,'c','Fontweight','bold');
grid on

