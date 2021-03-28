% Script to make quick map of dominant PFT and biomes (latter based on Smith et al., 2014, Biogeosciences, 11, 2027?2054)
%
% Dependencies:
% - readmasks_func.m
% - cbrewer (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
% - hickler_biome_read.m
%
%T. Pugh
%22.01.20

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

writetxt=false; %Write out to text file

obsbiome=true; %Whether to also plot observation-based biomes

outfile_name='simplemodel_best_est_100patch_10pCanopyCover_biomes_v2.txt';

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

addpath('/Users/pughtam/data/Hickler_vegmap')

%--- Read in data ---

lai=squeeze(lpj_to_grid_func_centre([lpjg_dir,'/lai_1961_1990'],1,0));
lai(:,:,13)=[]; %Remove total column

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);


%--- Calculate the dominant PFT by LAI ---
pftnames={'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G'};

dompft=NaN(360,720);
for ii=1:720
    for jj=1:360
        temp=lai(jj,ii,:);
        if max(temp)>0.0
            aa=find(temp==max(temp));
            if ~isempty(aa)
                if length(aa)>1
                    %HACK: Just take the first value (better to separate by NPP
                    %in this case)
                    dompft(jj,ii)=aa(1);
                else
                    dompft(jj,ii)=aa;
                end
            end
        end
    end
end
clear ii jj aa


%--- Calculate the biomes ---
% Calculate the biomes following the approach in Smith et al. (2014, Biogeosciences, 11, 2027?2054), editing code directly
% from the "biomes" function in the LPJ-GUESS code (v4.1).

%Some sums needed in the classification
treelai=nansum(lai(:,:,1:10),3);
grasslai=nansum(lai(:,:,11:12),3);
totlai=nansum(lai,3);
trbelai=nansum(lai(:,:,8:9),3);
bnelai=nansum(lai(:,:,1:2),3);
trlai=nansum(lai(:,:,8:10),3); %Tropical trees
telai=nansum(lai(:,:,4:7),3); %Temperate trees
btlai=nansum(lai(:,:,1:3),3); %Boreal trees
bneelai=nansum(lai(:,:,1:2),3);

lats=repmat(-89.75:0.5:89.75,[720 1])';

biomenames={'Boreal decid forest','Boreal ever forest','Temp/boreal mix fo.','Temp conifer forest','Temp decid forest',...
    'Temp broad ever fo.','Temp mixed forest','Trop season forest','Trop rain forest','Trop decid forest',...
    'Moist savannas','Dry savannas','Tall grassland','Dry grassland','Xeric wood/shrub','Arid shrub/steppe',...
    'Desert','Arctic/alpine tundra'};

pftnames={'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G'};

%NOTE: Below differs from Smith et al. (2014) in that treelai>2, rather than >2.5.
%treelaithres=2.0;
treelaithres=2.5;
biome=NaN(360,720);
for xx=1:720
    for yy=1:360
        if (treelai(yy,xx) > treelaithres) && (trbelai(yy,xx) > (0.6*treelai(yy,xx)))
            biome(yy,xx)=9; %Trop rain forest
        elseif (treelai(yy,xx) > treelaithres && (lai(yy,xx,10) > 0.6*treelai(yy,xx)))
            biome(yy,xx)=10; %Trop decid forest
        elseif (treelai(yy,xx) > treelaithres && (trlai(yy,xx) > 0.5*treelai(yy,xx)) && ...
                ((trbelai(yy,xx) > lai(yy,xx,7) && trbelai(yy,xx) > lai(yy,xx,5)) ||...
                (lai(yy,xx,10) > lai(yy,xx,7) && lai(yy,xx,10) > lai(yy,xx,5))))
            biome(yy,xx)=8; %Trop season forest
        elseif (treelai(yy,xx) > treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) &&...
                ((bneelai(yy,xx) > lai(yy,xx,3)) || (lai(yy,xx,6) > lai(yy,xx,3)))
            biome(yy,xx)=2; %Boreal ever forest
        elseif (treelai(yy,xx) > treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) &&...
                (lai(yy,xx,3) > bneelai(yy,xx)) && (lai(yy,xx,3) > lai(yy,xx,6))
            biome(yy,xx)=1; %Boreal decid forest
        elseif (treelai(yy,xx) > treelaithres) && (telai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,7) > 0.5*treelai(yy,xx))
            biome(yy,xx)=6; %Temp broad ever fo.
        elseif (treelai(yy,xx) > treelaithres) && (telai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,5) > 0.5*treelai(yy,xx))
            biome(yy,xx)=5; %Temp decid forest
        elseif (treelai(yy,xx) > treelaithres) && (telai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,4) > 0.5*treelai(yy,xx))
            biome(yy,xx)=4; %Temp conifer forest
        elseif (treelai(yy,xx) > treelaithres) && (btlai(yy,xx) > 0.2*treelai(yy,xx))
            biome(yy,xx)=3; %Temp/boreal mix fo.
        elseif (treelai(yy,xx) > treelaithres)
            biome(yy,xx)=7; %Temp mixed forest
        elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) && (bneelai(yy,xx) > lai(yy,xx,3) || lai(yy,xx,6) > lai(yy,xx,3))
            biome(yy,xx)=2; %Boreal ever forest
        elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,3) > bneelai(yy,xx) && lai(yy,xx,3) > lai(yy,xx,6))
            biome(yy,xx)=1; %Boreal decid forest
        elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (treelai(yy,xx) > 0.8*totlai(yy,xx))
            biome(yy,xx)=15; %Xeric wood/shrub
        elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (totlai(yy,xx) > 2.0)
            biome(yy,xx)=11; %Moist savannas
        elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres)
            biome(yy,xx)=12; %Dry savannas
        elseif (treelai(yy,xx) < 0.5) && (grasslai(yy,xx) > 0.2) && (lats(yy,xx) > 54)
            biome(yy,xx)=18; %Arctic/alpine tundra
        elseif (grasslai(yy,xx) > 2.0)
            biome(yy,xx)=13; %Tall grassland
        elseif (treelai(yy,xx) > 0.2) && (grasslai(yy,xx) < 1.0)
            biome(yy,xx)=16; %Arid shrub/steppe
        elseif (grasslai(yy,xx) > 0.2)
            biome(yy,xx)=14; %Dry grassland
        elseif (totlai(yy,xx) > 0.2)
            biome(yy,xx)=16; %Arid shrub/steppe
        elseif (totlai(yy,xx) <= 0.2)
            biome(yy,xx)=17; %Desert
        end
    end
end

%Read the observation based biomes?
if obsbiome
    [hickler_biomes,hickler_biomenames]=hickler_biome_read(false);
end

if use_cvegmask
    dompft=dompft.*cvegmask;
    biome=biome.*cvegmask;
    if obsbiome
        hickler_biomes=hickler_biomes.*cvegmask;
    end
end

if use_fmask
    dompft=dompft.*fmask;
    biome=biome.*fmask;
    if obsbiome
        hickler_biomes=hickler_biomes.*fmask;
    end
end

if use_bmask
    dompft=dompft.*bmask;
    biome=biome.*bmask;
    if obsbiome
        hickler_biomes=hickler_biomes.*bmask;
    end
end


%--- Make the map for dominant PFT ---

%Make figure
cmin=1;
cmax=12;

%Read mask for ocean areas
load(ocean_file);
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
oceanm=oceanm';

lons=-180:0.5:179.5;
lats=-90:0.5:89.5;

figure
%cmap=colormap(lines(13));
cmap=cbrewer('qual','Paired',13);
cmap=[0.9 0.9 0.9; cmap];

colormap(cmap)
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,dompft);
set(p1,'linestyle','none')
axis tight
caxis([-1 12])
c1=colorbar;
set(c1,'FontSize',12,'FontWeight','Bold')
set(c1,'Limits',[1 cmax])
set(c1,'Ticks',1:12,'TickLabels',pftnames)


%--- Make the map for biomes ---

%Make figure
cmin=1;
cmax=18;

%Read mask for ocean areas
load(ocean_file);
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
oceanm=oceanm';

lons=-180:0.5:179.5;
lats=-90:0.5:89.5;

figure
if obsbiome
    subplot(2,1,1)
end
cmap=cbrewer('qual','Paired',19);
cmap=[0.9 0.9 0.9; cmap];

colormap(cmap)
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,biome);
set(p1,'linestyle','none')
axis tight
caxis([-1 18])
c1=colorbar;
set(c1,'FontSize',12,'FontWeight','Bold')
set(c1,'Limits',[1 cmax])
set(c1,'Ticks',1:18,'TickLabels',biomenames)

if obsbiome
    subplot(2,1,2)
    
    hold on
    axesm('MapProjection','robinson','MapLatLimit',[23 80])
    hold on
    l1=pcolorm(lats,lons,oceanm);
    set(l1,'linestyle','none')
    p1=pcolorm(lats,lons,hickler_biomes);
    set(p1,'linestyle','none')
    axis tight
    caxis([-1 18])
    c1=colorbar;
    set(c1,'FontSize',12,'FontWeight','Bold')
    set(c1,'Limits',[1 cmax])
    set(c1,'Ticks',1:18,'TickLabels',hickler_biomenames)
end

%--- Write out to text file ---

if writetxt
    biome_out=biome;
    biome_out(isnan(biome))=-9999;
    biome_out=int32(biome_out);
    
    fid=fopen(outfile_name,'w');
    for yy=360:-1:1
        fprintf(fid,repmat('%7d',1,720),biome_out(yy,:));
        fprintf(fid,'\n');
    end
    clear yy
    fclose(fid);
end