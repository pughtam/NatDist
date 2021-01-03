% Plot biomes from LPJ-GUESS in comparison to those from observation-based estimates.
%
% T. Pugh
% 03.01.21

use_cvegmask=false; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

%writetxt=false; %Write out to text file
%outfile_name='simplemodel_best_est_100patch_10pCanopyCover_biomes_v2.txt';

%lpjg_dir='/Users/pughtam/LPJG/disturbance_prognostic_runs/simplemodel_best_est_100patch';
%lpjg_dir='/Users/pughtam/LPJG/bugfix_lu_progdist_ageout_runs/best_est';
%lpjg_dir='/Users/pughtam/LPJG/bugfix_lu_progdist_ageout_runs/best_est_luh2';
%lpjg_dir='/Users/pughtam/LPJG/bugfix_lu_progdist_ageout_runs/best_est_adjparam';
lpjg_dir='/Users/pughtam/LPJG/bugfix_lu_progdist_ageout_runs/best_est_adjparam_latosa4000';
%lpjg_dir='/Users/pughtam/LPJG/trunk_r8640_benchmark';
%lpjg_dir='/Users/pughtam/LPJG/trunk_r8278_aitkinresp/Output_global';
%lpjg_dir='/Users/pughtam/LPJG/disturbance_prognostic_runs/100_flat_orig_r5511/postproc/';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file

addpath('../obs_biomes/')


%--- Read in data ---

lai=squeeze(lpj_to_grid_func_centre([lpjg_dir,'/lai_1961_1990'],1,0));
lai(:,:,13)=[]; %Remove total column

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


%--- Calculate the dominant PFT by LAI ---

[biome,biomenames]=lpjg_biome_func(lai,2);

[biome_hengl,biomenames_hengl]=hengl_biome_read(true);

[biome_hickler,biomenames_hickler]=hickler_biome_read(2);

if use_cvegmask
    biome_hengl=biome_hengl.*cvegmask;
    biome_hickler=biome_hickler.*cvegmask;
    biome=biome.*cvegmask;
end

if use_fmask
    biome_hengl=biome_hengl.*fmask;
    biome_hickler=biome_hickler.*fmask;
    biome=biome.*fmask;
end

if use_bmask
    biome_hengl=biome_hengl.*bmask;
    biome_hickler=biome_hickler.*bmask;
    biome=biome.*bmask;
end

% Remove the tropical savannahs for plotting
biome_hengl_plot=biome_hengl;
biome_hengl_plot(biome_hengl_plot==10)=NaN;
biome_hengl_plot(biome_hengl_plot>10)=biome_hengl_plot(biome_hengl_plot>10)-1;
biome_hickler_plot=biome_hickler;
biome_hickler_plot(biome_hickler_plot==10)=NaN;
biome_hickler_plot(biome_hickler_plot>10)=biome_hickler_plot(biome_hickler_plot>10)-1;
biome_plot=biome;
biome_plot(biome_plot==10)=NaN;
biome_plot(biome_plot>10)=biome_plot(biome_plot>10)-1;

% Remove the tropical forest biomes for plotting
biome_hengl_plot=biome_hengl_plot-3;
biome_hengl_plot(biome_hengl_plot<1)=NaN;
biome_hickler_plot=biome_hickler_plot-3;
biome_hickler_plot(biome_hickler_plot<1)=NaN;
biome_plot=biome_plot-3;
biome_plot(biome_plot<1)=NaN;


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
s1=subplot(3,1,1);
cmap=cbrewer('qual','Paired',9);
cmap=[0.9 0.9 0.9; cmap];

colormap(cmap)
hold on
a1=axesm('MapProjection','robinson','MapLatLimit',[23 80]);
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,biome_plot);
set(p1,'linestyle','none')
axis tight
caxis([-1 9])
t1=textm(28,-175,'(a) LPJ-GUESS');
set(t1,'FontSize',12,'FontWeight','Bold')
% c1=colorbar;
% set(c1,'FontSize',12,'FontWeight','Bold')
% set(c1,'Limits',[1 cmax])
% set(c1,'Ticks',1.5:9,'TickLabels',biomenames([4:9,11:13]))

s2=subplot(3,1,2);
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l2=pcolorm(lats,lons,oceanm);
set(l2,'linestyle','none')
p2=pcolorm(lats,lons,biome_hengl_plot);
set(p2,'linestyle','none')
axis tight
caxis([-1 9])
t2=textm(28,-175,'(b) Hengl et al. (2018)');
set(t2,'FontSize',12,'FontWeight','Bold')
% c2=colorbar;
% set(c2,'FontSize',12,'FontWeight','Bold')
% set(c2,'Limits',[1 cmax])
% set(c2,'Ticks',1.5:9.5,'TickLabels',biomenames_hengl([4:9,11:13]))

s3=subplot(3,1,3);
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l3=pcolorm(lats,lons,oceanm);
set(l3,'linestyle','none')
p3=pcolorm(lats,lons,biome_hickler_plot);
set(p1,'linestyle','none')
axis tight
caxis([-1 9])
t3=textm(28,-175,'(c) Haxel. & Pren. (1996)');
set(t3,'FontSize',12,'FontWeight','Bold')
c3=colorbar;
set(c3,'FontSize',12,'FontWeight','Bold')
set(c3,'Limits',[1 cmax])
set(c3,'Ticks',1.5:9.5,'TickLabels',biomenames_hickler([4:9,11,13]))

set(s1,'Position',[0.01 0.67 0.7750 0.3])
set(s2,'Position',[0.01 0.34 0.7750 0.3])
set(s3,'Position',[0.01 0.01 0.7750 0.3])
set(c3,'Position',[0.785 0.34 0.0148 0.3])


%--- Write out to text file ---

% if writetxt
%     biome_out=biome;
%     biome_out(isnan(biome))=-9999;
%     biome_out=int32(biome_out);
%     
%     fid=fopen(outfile_name,'w');
%     for yy=360:-1:1
%         fprintf(fid,repmat('%7d',1,720),biome_out(yy,:));
%         fprintf(fid,'\n');
%     end
%     clear yy
%     fclose(fid);
% end

