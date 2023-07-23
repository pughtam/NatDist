function [biomes_0p5,biomenames]=hengl_biome_read(simplebiomes)
% Read Hengl et al. (2018) potential natural vegetation biome distribution.
% Aggregate to 0.5 x 0.5 degree resolution
% Optionally aggregate to groups for comparison with LPJ-GUESS forest biomes
%
% Dataset from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QQHCIK
%
% Original biome classification (LPJ-GUESS mapping in parentheses)
% 1     tropical evergreen broadleaf forest (tropical rain forest)
% 2     tropical semi-evergreen broadleaf forest (tropical seasonal forest)
% 3     tropical deciduous broadleaf forest and woodland (tropical deciduous forest)
% 4     warm-temperate evergreen and mixed forest (temperate broadleaved evergreen forest & temperate mixed forest)
% 7     cool-temperate rainforest (temperate/boreal mixed forest)
% 8     cool evergreen needleleaf forest (temperate/boreal mixed forest) (Temperate/boreal needleleaf evergreen forest)
% 9     cool mixed forest (temperate/boreal mixed forest)
% 13	temperate deciduous broadleaf forest (temperate deciduous forest)
% 14	cold deciduous forest (boreal deciduous forest) (boreal needleleaf deciduous forest)
% 15	cold evergreen needleleaf forest (boreal evergreen forest) (Temperate/boreal needleleaf evergreen forest)
% 16	temperate sclerophyll woodland and shrubland (xeric woodland/shrubland)
% 17	temperate evergreen needleleaf open woodland (xeric woodland/shrubland)
% 18	tropical savanna (dry savannah, moist savannah) (Tropical savannah) *USE A C4 CONDITION TO SEPARATE FROM STEPPE
% 19	<empty>
% 20	xerophytic woods/scrub (xeric woodland/shrubland)
% 22	steppe (tall grassland, dry grassland, arid woodland/steppe) (Grassland)
% 27	desert (desert)
% 28	graminoid and forb tundra (arctic/alpine tundra)
% 30	erect dwarf shrub tundra (arctic/alpine tundra)
% 31	low and high shrub tundra (arctic/alpine tundra)
% 32	prostrate dwarf shrub tundra (arctic/alpine tundra)
%
% T. Pugh
% 05.07.20

biomes=readgeoraster('/Users/pughtam/data/hengl_biomes/pnv_biome.type_biome00k_c_1km_s0..0cm_2000..2017_v0.1.tif');
biomes=flipud(biomes);

ndims=size(biomes);
cellsize=360/ndims(2);

% Optionally create simplified classification for comparison with LPJ-GUESS
if simplebiomes
    biomes_sim=biomes;
    %Numbers of 1-3 unchanged
    biomes_sim(biomes==4)=4; %Temperate evergreen forest
    biomes_sim(biomes==7 | biomes==9)=5; %Temperate/boreal mixed forest
    biomes_sim(biomes==8 | biomes==15)=6; %Temperate/boreal needleleaf evergreen forest
    biomes_sim(biomes==13)=7; %Temperate deciduous forest
    biomes_sim(biomes==14)=8; %Boreal deciduous forest
    biomes_sim(biomes==16 | biomes==17 | biomes==20)=9; %Xeric woodland/shrubland
    biomes_sim(biomes==18)=10; %Tropical savannah
    biomes_sim(biomes==22)=11; %Grassland
    biomes_sim(biomes==27)=12; %Desert
    biomes_sim(biomes==28 | biomes==30 | biomes==31 | biomes==32)=13; %Tundra
    
    biomes=biomes_sim;
    clear biomes_sim
    
    biomenames={'Tropical evergreen forest','Tropical seasonal forest','Tropical deciduous forest',...
        'Temperate evergreen forest','Temperate/boreal mixed forest','Needleleaf evergreen forest',...
        'Temperate deciduous forest','Boreal deciduous forest',...
        'Xeric woodland/shrubland','Tropical savannah','Grassland','Desert','Tundra'};
else
    biomenames={'tropical evergreen broadleaf forest','tropical semi-evergreen broadleaf forest',...
        'tropical deciduous broadleaf forest and woodland,','warm-temperate evergreen and mixed forest',...
        'cool-temperate rainforest','cool evergreen needleleaf forest','cool mixed forest',...
        'temperate deciduous broadleaf forest','cold deciduous forest','cold evergreen needleleaf forest',...
        'temperate sclerophyll woodland and shrubland','temperate evergreen needleleaf open woodland',...
        'tropical savanna','xerophytic woods/scrub','steppe','desert','graminoid and forb tundra',...
        'erect dwarf shrub tundra','low and high shrub tundra','prostrate dwarf shrub tundra'};
end

% Aggregate to 0.5 degrees
ncells=0.5/cellsize;
biomes_0p5=NaN(360,720);
for xx=1:720
    for yy=56:353
        yyy=yy-55;
        xx_s=(xx*ncells)-ncells+1;
        xx_e=xx*ncells;
        yy_s=(yyy*ncells)-ncells+1;
        yy_e=yyy*ncells;
        temp=biomes(yy_s:yy_e,xx_s:xx_e);
        temp(temp==255)=NaN;
        biomes_0p5(yy,xx)=mode(temp(:));
    end
end
clear xx yy yyy xx_s xx_e yy_s yy_e temp
biomes_0p5(biomes_0p5==0)=NaN;

if ~simplebiomes
    %Re-number to continuous values
    usednums=unique(biomes_0p5(isnan(biomes_0p5)==0));
    biomes_0p5_renum=NaN(size(biomes_0p5));
    for nn=1:length(usednums)
        biomes_0p5_renum(biomes_0p5==usednums(nn))=nn;
    end
    clear nn
    biomes_0p5=biomes_0p5_renum;
    clear biomes_0p5_renum
end
