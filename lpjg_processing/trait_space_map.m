% Make a map of all grid cells that in the trait space of the 77 landscapes
%
% T. Pugh (based on R function by C. Senf)
% 10.08.23

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
use_bmask=true; %Mask by temperate/boreal biomes

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';

fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='./data/';
ocean_file='/Users/pughtam/data/ESA_landcover/esa_05_landcover.mat'; %Ocean mask file

% Get the LPJ-GUESS simulated climate and trait data for the disturbance model
temprange=permute(ncread([lpjg_dir,'/temprange_LPJ-GUESS_standard_nat_2014.nc'],'temprange'),[2 1]); 
temprange(temprange==0)=NaN;

wooddensity=permute(ncread([lpjg_dir,'/wooddensity_LPJ-GUESS_standard_nat_2014.nc'],'wooddensity'),[2 1]); 
wooddensity(wooddensity==0)=NaN;

% Get the landscape-level climate and disturbance data
data_landscape=readtable('../disturbance_rates/data/traits_landscape.csv');

% Define the trait space occupied
chull=convhull(data_landscape.temp_range,data_landscape.wooddensity_species);

% Make a check plot for the convex hull
figure
plot(data_landscape.temp_range,data_landscape.wooddensity_species,'.')
hold on
plot(data_landscape.temp_range(chull),data_landscape.wooddensity_species(chull))

inrange=false(size(temprange));
for xx=1:720
    for yy=1:360
        if ~isnan(temprange(yy,xx))
            if inpolygon(temprange(yy,xx),wooddensity(yy,xx),...
                    data_landscape.temp_range(chull),data_landscape.wooddensity_species(chull))
                inrange(yy,xx)=true;
            end
        end
    end
end
clear yy xx

% Apply masks
[cvegmask,fmask,bmask,ffrac]=readmasks_func(use_cvegmask,use_fmask,false,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir);

if use_cvegmask
    inrange=inrange.*cvegmask;
end

if use_fmask
    inrange=inrange.*fmask;
end

if use_bmask
    inrange=inrange.*bmask;
end

% Make a map of the grid cells which are in the rnage of the landscape trait data.
load(ocean_file);
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
oceanm=oceanm';

lons=-180:0.5:179.5;
lats=-90:0.5:89.5;

figure
cmap=[0.9 0.9 0.9; 0.5 0 0; 0 0.5 0];

colormap(cmap)
axesm('MapProjection','robinson','MapLatLimit',[23 80])
hold on
l1=pcolorm(lats,lons,oceanm);
set(l1,'linestyle','none')
p1=pcolorm(lats,lons,inrange);
set(p1,'linestyle','none')
axis tight

