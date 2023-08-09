% Plot biomes from LPJ-GUESS in comparison to those from observation-based estimates.
%
% Dependencies
% - lpjg_biome_func-m
% - hengl_biome_read.m
% - hickler_biome_read.m
% - readmasks_func.m
% - lpj_to_grid_func_centre.m
% - cbrewer (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
%
% T. Pugh
% 03.01.21

use_cvegmask=false; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='./data/';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file

addpath('../obs_biomes/')


%--- Read in data ---

lai_1930=permute(ncread([lpjg_dir,'/LAI_LPJ-GUESS_standard_nat_1930.nc'],'LAI'),[2 1 3]);

lai_2014=permute(ncread([lpjg_dir,'/LAI_LPJ-GUESS_standard_nat_2014.nc'],'LAI'),[2 1 3]);

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


%--- Calculate the dominant PFT by LAI ---

[biome_1930,biomenames]=lpjg_biome_func(lai_1930,2);
[biome_2014,biomenames]=lpjg_biome_func(lai_2014,2);

[biome_hengl,biomenames_hengl]=hengl_biome_read(true);

[biome_hickler,biomenames_hickler]=hickler_biome_read(2);

if use_cvegmask
    biome_hengl=biome_hengl.*cvegmask;
    biome_hickler=biome_hickler.*cvegmask;
    biome_1930=biome_1930.*cvegmask;
    biome_2014=biome_2014.*cvegmask;
end

if use_fmask
    biome_hengl=biome_hengl.*fmask;
    biome_hickler=biome_hickler.*fmask;
    biome_1930=biome_1930.*fmask;
    biome_2014=biome_2014.*fmask;
end

if use_bmask
    biome_hengl=biome_hengl.*bmask;
    biome_hickler=biome_hickler.*bmask;
    biome_1930=biome_1930.*bmask;
    biome_2014=biome_2014.*bmask;
end

% Remove the tropical savannahs for plotting
biome_hengl_plot=biome_hengl;
biome_hengl_plot(biome_hengl_plot==10)=NaN;
biome_hengl_plot(biome_hengl_plot>10)=biome_hengl_plot(biome_hengl_plot>10)-1;
biome_hickler_plot=biome_hickler;
biome_hickler_plot(biome_hickler_plot==10)=NaN;
biome_hickler_plot(biome_hickler_plot>10)=biome_hickler_plot(biome_hickler_plot>10)-1;
biome_1930_plot=biome_1930;
biome_1930_plot(biome_1930_plot==10)=NaN;
biome_1930_plot(biome_1930_plot>10)=biome_1930_plot(biome_1930_plot>10)-1;
biome_2014_plot=biome_2014;
biome_2014_plot(biome_2014_plot==10)=NaN;
biome_2014_plot(biome_2014_plot>10)=biome_2014_plot(biome_2014_plot>10)-1;

% Remove the tropical forest biomes for plotting
biome_hengl_plot=biome_hengl_plot-3;
biome_hengl_plot(biome_hengl_plot<1)=NaN;
biome_hickler_plot=biome_hickler_plot-3;
biome_hickler_plot(biome_hickler_plot<1)=NaN;
biome_1930_plot=biome_1930_plot-3;
biome_1930_plot(biome_1930_plot<1)=NaN;
biome_2014_plot=biome_2014_plot-3;
biome_2014_plot(biome_2014_plot<1)=NaN;


%--- Make the map for biomes ---

%Make figure
cmin=1;
cmax=9; %13;

%Read mask for ocean areas
load(ocean_file);
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
oceanm=oceanm';

lons=-180:0.5:179.5;
lats=-90:0.5:89.5;

figure
s1=subplot(4,1,1);
cmap=cbrewer('qual','Paired',9);
cmap=[0.9 0.9 0.9; cmap];

colormap(cmap)
hold on
a1=axesm('MapProjection','robinson','MapLatLimit',[23 80]);
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,biome_2014_plot);
set(p1,'linestyle','none')
axis tight
caxis([-1 9])
t1=textm(28,-175,'(a) LPJ-GUESS, 2001-2014');
set(t1,'FontSize',12,'FontWeight','Bold')

s2=subplot(4,1,2);
hold on
a2=axesm('MapProjection','robinson','MapLatLimit',[23 80]);
hold on
l2=pcolorm(lats,lons,oceanm);
set(l2,'linestyle','none')
p2=pcolorm(lats,lons,biome_1930_plot);
set(p2,'linestyle','none')
axis tight
caxis([-1 9])
t2=textm(28,-175,'(b) LPJ-GUESS, 1901-1930');
set(t2,'FontSize',12,'FontWeight','Bold')

s3=subplot(4,1,3);
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l3=pcolorm(lats,lons,oceanm);
set(l3,'linestyle','none')
p3=pcolorm(lats,lons,biome_hengl_plot);
set(p3,'linestyle','none')
axis tight
caxis([-1 9])
t3=textm(28,-175,'(c) Hengl et al. (2018)');
set(t3,'FontSize',12,'FontWeight','Bold')

s4=subplot(4,1,4);
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l4=pcolorm(lats,lons,oceanm);
set(l4,'linestyle','none')
p4=pcolorm(lats,lons,biome_hickler_plot);
set(p4,'linestyle','none')
axis tight
caxis([-1 9])
t4=textm(28,-175,'(d) Haxel. & Pren. (1996)');
set(t4,'FontSize',12,'FontWeight','Bold')
c4=colorbar;
set(c4,'FontSize',12,'FontWeight','Bold')
set(c4,'Limits',[1 cmax])
set(c4,'Ticks',1.5:9.5,'TickLabels',biomenames_hickler([4:9,11,13]))

set(s1,'Position',[0.01 0.75 0.7750 0.23])
set(s2,'Position',[0.01 0.50 0.7750 0.23])
set(s3,'Position',[0.01 0.25 0.7750 0.23])
set(s4,'Position',[0.01 0.01 0.7750 0.23])
set(c4,'Position',[0.8 0.34 0.0148 0.3])

