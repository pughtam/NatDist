function [biome,biomenames]=lpjg_biome_func(lai,biome_scheme)
% Function to convert LPJ-GUESS LAI output into a biome classification
% - 1. Smith et al. (2014, Biogeosciences, 11, 2027-2054)
% - 2. Hengl et al. (2018) biomes, as preprocessed in hengl_biome_read.m
% Note that the settings herein assume an order of PFTs in the lai array passed to this function as follows:
% 'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G'
%
% T. Pugh
% 05.07.20


%Some sums needed in the classification
treelai=nansum(lai(:,:,1:10),3);
grasslai=nansum(lai(:,:,11:12),3);
totlai=nansum(lai,3);
trbelai=nansum(lai(:,:,8:9),3); %TrBE & TrIBE
trlai=nansum(lai(:,:,8:10),3); %Tropical trees
telai=nansum(lai(:,:,4:7),3); %Temperate trees
btlai=nansum(lai(:,:,1:3),3); %Boreal trees
bneelai=nansum(lai(:,:,1:2),3); %BNE & BINE

lats=repmat(-89.75:0.5:89.75,[720 1])';

if biome_scheme==1
    % Calculate the biomes following the approach in Smith et al. (2014, Biogeosciences, 11, 2027-2054), editing code directly
    % from the "biomes" function in the LPJ-GUESS code (v4.1).
    
    biomenames={'Boreal decid forest','Boreal ever forest','Temp/boreal mix fo.','Temp conifer forest','Temp decid forest',...
        'Temp broad ever fo.','Temp mixed forest','Trop season forest','Trop rain forest','Trop decid forest',...
        'Moist savannas','Dry savannas','Tall grassland','Dry grassland','Xeric wood/shrub','Arid shrub/steppe',...
        'Desert','Arctic/alpine tundra'};
    
    %treelaithres=2.0; %NOTE: differs from Smith et al. (2014) in that treelai>2, rather than >2.5.
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
    
    biomenames={'Trop rain forest','Trop seasonal forest','Trop decid forest',...
        'Temp ever fo.','Temp/boreal mix fo.','Temp/Boreal ever forest','Temp decid forest',...
        'Boreal decid forest','Xeric wood/shrub','Tropical savannas','Grassland',...
        'Desert','Tundra'};
    
    %NOTE: Below differs from Smith et al. (2014) in that treelai>2, rather than >2.5.
    treelaithres=2.0;
    biome=NaN(360,720);
    for xx=1:720
        for yy=1:360
            if (treelai(yy,xx) > treelaithres) && (trbelai(yy,xx) > (0.6*treelai(yy,xx)))
                biome(yy,xx)=1; %7; %Trop rain forest
            elseif (treelai(yy,xx) > treelaithres && (lai(yy,xx,10) > 0.6*treelai(yy,xx)))
                biome(yy,xx)=3; %8; %Trop decid forest
            elseif (treelai(yy,xx) > treelaithres && (trlai(yy,xx) > 0.5*treelai(yy,xx)) && ...
                    ((trbelai(yy,xx) > lai(yy,xx,7) && trbelai(yy,xx) > lai(yy,xx,5)) ||...
                    (lai(yy,xx,10) > lai(yy,xx,7) && lai(yy,xx,10) > lai(yy,xx,5))))
                biome(yy,xx)=2; %6; %Trop season forest
            elseif (treelai(yy,xx) > treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) &&...
                    ((bneelai(yy,xx) > lai(yy,xx,3)) || (lai(yy,xx,6) > lai(yy,xx,3)))
                biome(yy,xx)=6; %2; %Temperate/boreal needleleaf evergreen forest
            elseif (treelai(yy,xx) > treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) &&...
                    (lai(yy,xx,3) > bneelai(yy,xx)) && (lai(yy,xx,3) > lai(yy,xx,6))
                biome(yy,xx)=8; %1; %Boreal deciduous forest
            elseif (treelai(yy,xx) > treelaithres) && (telai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,7) > 0.5*treelai(yy,xx))
                biome(yy,xx)=4; %5; %Temperate evergreen forest
            elseif (treelai(yy,xx) > treelaithres) && (telai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,5) > 0.5*treelai(yy,xx))
                biome(yy,xx)=7; %4; %Temperate deciduous forest
            elseif (treelai(yy,xx) > treelaithres) && (telai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,4) > 0.5*treelai(yy,xx))
                biome(yy,xx)=6; %2; %Temperate/boreal needleleaf evergreen forest
            elseif (treelai(yy,xx) > treelaithres)
                biome(yy,xx)=5; %3; %Temperate/boreal mixed forest 
            elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) && (bneelai(yy,xx) > lai(yy,xx,3) || lai(yy,xx,6) > lai(yy,xx,3))
                biome(yy,xx)=6; %2; %Temperate/boreal needleleaf evergreen forest
            elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (btlai(yy,xx) > 0.8*treelai(yy,xx)) && (lai(yy,xx,3) > bneelai(yy,xx) && lai(yy,xx,3) > lai(yy,xx,6))
                biome(yy,xx)=8; %1; %Boreal deciduous forest
            elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (treelai(yy,xx) > 0.8*totlai(yy,xx))
                biome(yy,xx)=9; %11; %Xeric wood/shrub
            elseif (treelai(yy,xx) > 0.5) && (treelai(yy,xx) < treelaithres) && (lai(yy,xx,12) > 0.8*grasslai(yy,xx))
                biome(yy,xx)=10; %9; %Tropical savannah
            elseif (treelai(yy,xx) < 0.5) && (grasslai(yy,xx) > 0.2) && (lats(yy,xx) > 54)
                biome(yy,xx)=13; %Tundra
            elseif (totlai(yy,xx) > 2.0)
                biome(yy,xx)=11; %10; %Grassland
            elseif (totlai(yy,xx) <= 0.2)
                biome(yy,xx)=12; %Desert
            end
        end
    end
    
else
    error('biome_scheme must be set to 1 (Smith) or 2 (Hengl)')
end
