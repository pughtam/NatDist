function [out_map,names]=hickler_biome_read(biome_scheme)

% 1 Boreal deciduous forest/woodland
% 2 Boreal evergreen forest/woodland
% 3 Temperate/boreal mixed forest
% 4 Temperate conifer forest
% 5 Temperate deciduous forest
% 6 Temperate broad-leaved evergreen forest
% 7 Temperate mixed forest (missing from documentation, but consistent with papers)
% 8 Tropical seasonal forest
% 9 Tropical rain forest
% 10 Tropical deciduous Forest
% 11 Moist savannas
% 12 Dry savannas    
% 13 Tall grassland 
% 14 Short grassland
% 15 Xeric woodlands/scrub
% 16 Arid shrubland/steppe
% 17 Desert 
% 18 Arctic/alpine tundra

data_in=dlmread('vegmap18.out','',1,0);
lon_in=data_in(:,1);
lat_in=data_in(:,2);
bio_in=data_in(:,3);

lon_map=-179.75:0.5:179.75;
lat_map=-89.75:0.5:89.75;

bio_map=NaN(length(lat_map),length(lon_map));
for i=1:length(lon_in)
    lon_ind=find(lon_map-0.25 == (lon_in(i)/10));
    lat_ind=find(lat_map-0.25 == (lat_in(i)/10));
    
    bio_map(lat_ind,lon_ind)=bio_in(i);
end

if biome_scheme==1
    names={'Boreal deciduous forest/woodland','Boreal evergreen forest/woodland','Temperate/boreal mixed forest','Temperate conifer forest','Temperate deciduous forest',...
        'Temperate broad-leaved evergreen forest','Temperate mixed forest','Tropical seasonal forest','Tropical rain forest','Tropical deciduous Forest','Moist savannas','Dry savannas',...
        'Tall grassland','Short grassland','Xeric woodlands/scrub','Arid shrubland/steppe','Desert','Arctic/alpine tundra','Polar desert'};
    
    out_map=bio_map;
elseif biome_scheme==2
    %New biome classification for consistency with Hengl et al. (2018) biomes (simplified version), as preprocessed in hengl_biome_read.m
    % 1     Tropical evergreen broadleaf forest (tropical rain forest)
    % 2     Tropical semi-evergreen broadleaf forest (tropical seasonal forest)
    % 3     Tropical deciduous broadleaf forest and woodland (tropical deciduous forest)
    % 4     Temperate evergreen forest
    % 5     Temperate/boreal mixed forest
    % 6     Temperate/boreal needleleaf evergreen forest
    % 7     Temperate deciduous forest
    % 8     Boreal deciduous forest
    % 9     Xeric woodland/shrubland
    % 10    Tropical savannah
    % 11    Grassland
    % 12    Desert
    % 13    Tundra
    
    bio_map_cond=NaN(size(bio_map));
    bio_map_cond(bio_map==9)=1; %Tropical evergreen broadleaf forest (tropical rain forest)
    bio_map_cond(bio_map==8)=2; %Tropical semi-evergreen broadleaf forest (tropical seasonal forest)
    bio_map_cond(bio_map==10)=3; %Tropical deciduous broadleaf forest and woodland (tropical deciduous forest)
    bio_map_cond(bio_map==6)=4; %Temperate evergreen forest
    bio_map_cond(bio_map==7 | bio_map==3)=5; %Temperate/boreal mixed forest
    bio_map_cond(bio_map==2 | bio_map==4)=6; %Temperate/boreal needleleaf evergreen forest
    bio_map_cond(bio_map==5)=7; %Temperate deciduous forest
    bio_map_cond(bio_map==1)=8; %Boreal deciduous forest/woodland
    bio_map_cond(bio_map==15 | bio_map==16)=9; %Xeric woodland/shrubland
    bio_map_cond(bio_map==11 | bio_map==12)=10; %Savannahs (Tropical savannah)
    bio_map_cond(bio_map==13 | bio_map==14)=11; %Grassland
    bio_map_cond(bio_map==17)=12; %Desert
    bio_map_cond(bio_map==18)=13; %Tundra
    
    names={'Tropical evergreen forest','Tropical seasonal forest','Tropical deciduous forest',...
        'Temperate evergreen forest','Temperate/boreal mixed forest','Needleleaf evergreen forest',...
        'Temperate deciduous forest','Boreal deciduous forest',...
        'Xeric woodland/shrubland','Savannah','Grassland','Desert','Tundra'};
    
    out_map=bio_map_cond;
end