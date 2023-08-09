% Script to make quick map of fraction of tree cover which is broadleaf for LPJ-GUESS and ESA data
%
% Dependencies:
% - readmasks_func.m
% - cbrewer (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
% - esa_broadleaf_frac_0p5deg.mat (from esa_broadleaf_frac_0p5deg.m)
%
%T. Pugh
%03.04.20

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

writetxt=false; %Write out to text file

outfile_name='best_est_adjparam_10pCanopyCover_broadleaf_frac.txt';
outfile_esa_name='esa_10pCanopyCover_broadleaf_frac.txt';

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';
esa_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/NatDist_working/ESA_processing';
fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='./data/';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file


%--- Read in data ---

lai=permute(ncread([lpjg_dir,'/LAI_LPJ-GUESS_standard_nat_2014.nc'],'LAI'),[2 1 3]); 

[cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);

load([esa_dir,'/esa_broadleaf_frac_0p5deg.mat']);

%--- Calculate the LAI of each tree type ---

lai_needleleaf=sum(lai(:,:,1:4),3);
lai_broadleaf=sum(lai(:,:,5:10),3);
broadleaf_frac=lai_broadleaf./(lai_broadleaf+lai_needleleaf);

if use_cvegmask
    broadleaf_frac=broadleaf_frac.*cvegmask;
    broadleaf_frac_esa=broadleaf_frac_esa.*cvegmask;
end

if use_fmask
    broadleaf_frac=broadleaf_frac.*fmask;
    broadleaf_frac_esa=broadleaf_frac_esa.*fmask;
end

if use_bmask
    broadleaf_frac=broadleaf_frac.*bmask;
    broadleaf_frac_esa=broadleaf_frac_esa.*bmask;
end

broadleaf_frac=broadleaf_frac*100; %Convert to percentage


%--- Make the maps ---

%Make figure
cmin=0;
cmax=100;
ncgrad=100; %Number of colour gradations

%Read mask for ocean areas
load(ocean_file);
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
oceanm=oceanm';

lons=-180:0.5:179.5;
lats=-90:0.5:89.5;

%LPJ-GUESS
figure
s1=subplot(2,1,1);
cmap=colormap(parula(ncgrad));
cmap=[0.9 0.9 0.9; cmap];
cminraw=cmin-(cmax/ncgrad);

colormap(cmap)
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,broadleaf_frac);
set(p1,'linestyle','none')
axis tight
caxis([cminraw cmax])
c1=colorbar;
set(c1,'FontSize',12,'FontWeight','Bold')
set(c1,'Limits',[cmin cmax])
title('(a) LPJ-GUESS')

%ESA
s2=subplot(2,1,2);
cmap=colormap(parula(ncgrad));
cmap=[0.9 0.9 0.9; cmap];
cminraw=cmin-(cmax/ncgrad);

colormap(cmap)
hold on
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,broadleaf_frac_esa);
set(p1,'linestyle','none')
axis tight
caxis([cminraw cmax])
title('(b) ESA')

set(s1,'Position',[0.01 0.50 0.7750 0.46])
set(s2,'Position',[0.01 0.01 0.7750 0.46])
set(c1,'Position',[0.8 0.34 0.0148 0.3])


%--- Write out to text files ---

if writetxt
    %LPJ-GUESS
    broadleaf_frac_out=broadleaf_frac;
    broadleaf_frac_out(isnan(broadleaf_frac))=-9999;
    broadleaf_frac_out=int32(broadleaf_frac_out);
    
    fid=fopen(outfile_name,'w');
    for yy=360:-1:1
        fprintf(fid,repmat('%7d',1,720),broadleaf_frac_out(yy,:));
        fprintf(fid,'\n');
    end
    clear yy
    fclose(fid);
    
    %ESA
    broadleaf_frac_esa_out=broadleaf_frac_esa;
    broadleaf_frac_esa_out(isnan(broadleaf_frac_esa))=-9999;
    broadleaf_frac_esa_out=int32(broadleaf_frac_esa_out);
    
    fid=fopen(outfile_esa_name,'w');
    for yy=360:-1:1
        fprintf(fid,repmat('%7d',1,720),broadleaf_frac_esa_out(yy,:));
        fprintf(fid,'\n');
    end
    clear yy
    fclose(fid);
end