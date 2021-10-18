%Script to read ESA landcover mask from http://maps.elie.ucl.ac.be/CCI/viewer/download.php
%and make masks for different landcover types.
%T. Pugh
%25.03.17

%ESA landcover codes:
% 0;No data;0;0;0
% 10;Cropland, rainfed;255;255;100
% 11;Herbaceous cover;255;255;100
% 12;Tree or shrub cover;255;255;0
% 20;Cropland, irrigated or post-flooding;170;240;240
% 30;Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%);220;240;100
% 40;Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%) ;200;200;100
% 50;Tree cover, broadleaved, evergreen, closed to open (>15%);0;100;0
% 60;Tree cover, broadleaved, deciduous, closed to open (>15%);0;160;0
% 61;Tree cover, broadleaved, deciduous, closed (>40%);0;160;0
% 62;Tree cover, broadleaved, deciduous, open (15-40%);170;200;0
% 70;Tree cover, needleleaved, evergreen, closed to open (>15%);0;60;0
% 71;Tree cover, needleleaved, evergreen, closed (>40%);0;60;0
% 72;Tree cover, needleleaved, evergreen, open (15-40%);0;80;0
% 80;Tree cover, needleleaved, deciduous, closed to open (>15%);40;80;0
% 81;Tree cover, needleleaved, deciduous, closed (>40%);40;80;0
% 82;Tree cover, needleleaved, deciduous, open (15-40%);40;100;0
% 90;Tree cover, mixed leaf type (broadleaved and needleleaved);120;130;0
% 100;Mosaic tree and shrub (>50%) / herbaceous cover (<50%);140;160;0
% 110;Mosaic herbaceous cover (>50%) / tree and shrub (<50%);190;150;0
% 120;Shrubland;150;100;0
% 121;Shrubland evergreen;120;75;0
% 122;Shrubland deciduous;150;100;0
% 130;Grassland;255;180;50
% 140;Lichens and mosses;255;220;210
% 150;Sparse vegetation (tree, shrub, herbaceous cover) (<15%);255;235;175
% 152;Sparse shrub (<15%);255;210;120
% 153;Sparse herbaceous cover (<15%);255;235;175
% 160;Tree cover, flooded, fresh or brakish water;0;120;90
% 170;Tree cover, flooded, saline water;0;150;120
% 180;Shrub or herbaceous cover, flooded, fresh/saline/brakish water;0;220;130
% 190;Urban areas;195;20;0
% 200;Bare areas;255;245;215
% 201;Consolidated bare areas;220;220;220
% 202;Unconsolidated bare areas;255;245;215
% 210;Water bodies;0;70;200
% 220;Permanent snow and ice;255;255;255

    esa=geotiffread('ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif');
    
    %Aggregate to 0.5 degrees, taking the most common landcover type in each gridcell
    %360 units per degree longitude and latitude
    
    nlatcell=180;
    nloncell=180;
    esa_05=NaN(360,720);
    for ii=1:720
        for jj=1:360
            indlat_s=(jj*nlatcell)-nlatcell+1;
            indlat_e=jj*nlatcell;
            indlon_s=(ii*nloncell)-nloncell+1;
            indlon_e=ii*nloncell;
            temp=esa(indlat_s:indlat_e,indlon_s:indlon_e);
            esa_05(jj,ii)=mode(temp(:));
            clear temp
        end
        fprintf('%d\n',ii)
    end
    esa_05=flipud(esa_05);
    
    save esa_05_landcover.mat esa_05
