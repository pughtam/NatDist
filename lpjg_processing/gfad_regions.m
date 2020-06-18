function [new_region,regionnames,regionnames_short]=gfad_regions(gfad_filepath_stan,onedeg)
%Script to take the GFAD regions and process them into something more
%accessible
%
%T. Pugh
%16.08.19

%gfad_filepath_stan='/Users/pughtam/data/GFAD_V1-1/GFAD_V1-1.nc';
gfad_region=ncread(gfad_filepath_stan,'adminunit');
gfad_lat=ncread(gfad_filepath_stan,'lat');
gfad_lon=ncread(gfad_filepath_stan,'lon');
[lats,lons]=meshgrid(gfad_lat,gfad_lon);

%Region indices from GFAD
tropics=1;
russia=2;
us=3:54;
alaska=20;
canada=55:69;
china=70:101;
europe=102:145;
japan=148;
korea=149:150;
new_zealand=151;
mongolia=152;
kazak=153;

regionnames={'Rest of world','West Russia','Mid Russia','East Russia','North-East US',...
    'South-East US','West US','Alaska','East Canada','West Canada','South/Mid. China',...
    'North China','South/Mid. Europe','North Europe','Japan and Korea','New Zealand'};

regionnames_short={'RestOfWorld','WestRussia','MidRussia','EastRussia','NorthEastUS',...
    'SouthEastUS','WestUS','Alaska','EastCanada','WestCanada','SouthMidChina',...
    'NorthChina','SouthMidEurope','NorthEurope','JapanKorea','NewZealand'};

new_region=zeros(size(gfad_region));
new_region(gfad_region==tropics)=1; %Rest of world
new_region(gfad_region==russia & lons<60 & lons>0)=2; %West Russia
new_region(gfad_region==russia & lons>=60 & lons <=120)=3; %Mid Russia
new_region(gfad_region==russia & (lons>120 | lons<0))=4; %East Russia
for nn=1:length(us)
    new_region(gfad_region==us(nn) & lons>-100 & lats>37)=5; %North-East US
    new_region(gfad_region==us(nn) & lons>-100 & lats<=37)=6; %South-East US
    new_region(gfad_region==us(nn) & lons<=-100)=7; %West US
end
new_region(gfad_region==alaska)=8;
for nn=1:length(canada)
    new_region(gfad_region==canada(nn) & lons>-100)=9; %East Canada
    new_region(gfad_region==canada(nn) & lons<=-100)=10; %West Canada
end
for nn=1:length(china)
    new_region(gfad_region==china(nn) & lats<40)=11; %South and Mid China
    new_region(gfad_region==china(nn) & lats>=40)=12; %North China
end
for nn=1:length(europe)
    new_region(gfad_region==europe(nn) & lats<55)=13; %South and Mid Europe
    new_region(gfad_region==europe(nn) & lats>=55)=14; %North Europe
end
new_region(gfad_region==japan)=15; %Merge Japan and Korea
for nn=1:length(korea)
    new_region(gfad_region==korea(nn))=15; %Merge Japan and Korea
end
new_region(gfad_region==new_zealand)=16;
new_region(gfad_region==mongolia)=3; %Merge Mongolia with Mid Russia (because of minimal data)
new_region(gfad_region==kazak)=3; %Merge Kazakstan with Mid Russia (because of minimal data)
new_region(gfad_region==tropics & lats<23 & lats>-23 & lons<-32)=17; %South America tropics
new_region(gfad_region==tropics & lats<23 & lats>-23 & lons>-32 & lons<55)=18; %Africa tropics
new_region(gfad_region==tropics & lats<23 & lats>-23 & lons>55)=19; %South-east Asia tropics
clear nn

%p1=pcolor(flipud(new_region'));
%set(p1,'linestyle','none')

new_region_1deg=NaN(360,180);
if onedeg
    %Aggregate to 1 degree gridcells
    for xx=1:360
        for yy=1:180
            ii_s=(xx*2)-1;
            ii_e=xx*2;
            jj_s=(yy*2)-1;
            jj_e=yy*2;
            temp=new_region(ii_s:ii_e,jj_s:jj_e);
            new_region_1deg(xx,yy)=mode(temp(:));
            clear temp ii_s ii_e jj_s jj_e
        end
        clear yy
    end
    clear xx
    new_region=new_region_1deg;
    clear new_region_1deg
end
            
