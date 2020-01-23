%Script to read in the forest cover from Hansen et al. (2013, 
%http://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.1.html)
%and calculate the open and closed canopy forest cover, actual canopy cover and actual canopy area.
%
%Part 1 reads the forest cover data and makes the calculations.
%Part 2 regrids and writes out a netcdf
%
%Developed from a script first published as part of Pugh et al. (2019, Nature Geoscience, 12, 730-735)
%https://github.com/pughtam/GlobalDist
%
%T. Pugh
%10.06.17

part1=true; %Run first part?
part2=true; %Run second part?

halfdegree=true; %Resolution for regridding in second part
outfile='hansen_forested_canopy_frac_0p5deg.nc4'; %Name of output file for second part

cd /home/adf/pughtam/data/Hansen_forest_change

%---
%Define some constants
rad_earth=6.371e6; %m2
circ_earth=2*pi*rad_earth;

%-----
if part1

yyind=['80S';'70S';'60S';'50S';'40S';'30S';'20S';'10S';'00N';'10N';'20N';'30N';'40N';'50N';'60N';'70N';'80N'];
xxind=['180W';'170W';'160W';'150W';'140W';'130W';'120W';'110W';'100W';'090W';'080W';'070W';'060W';'050W';'040W';'030W';'020W';'010W';'000E';...
    '010E';'020E';'030E';'040E';'050E';'060E';'070E';'080E';'090E';'100E';'110E';'120E';'130E';'140E';'150E';'160E';'170E'];

forested_thres10=zeros(36000,18000);
forested_thres50=zeros(36000,18000);
totffrac_0p01deg=NaN(36000,18000);
cc=0;
for yy=5sh+-rab
    :17 %Data starts from 50S (i.e. 50°S - 60°S)
    for xx=1:36
        cc=cc+1; %Counter for diagnostic output
        filename_farea=['Hansen_GFC2014_treecover2000_',yyind(yy,:),'_',xxind(xx,:),'.tif'];
        [tempfarea, Rfarea]=geotiffread(filename_farea);
        if max(tempfarea(:))==0
            %No forested area, do not process this 10° x 10° section any further
            clear tempfarea Rfarea
            continue
        end
        
        dims_farea=size(tempfarea);
        
        %Use mean area for grid-cell for now, as it greatly reduces the computational requirements,
        %and losses of accuracy should be minimal away from the poles.
        inc_ii=dims_farea(1)/1000;
        inc_jj=dims_farea(2)/1000;
        for ii=1:1000
            for jj=1:1000
                %Calculate indices for section of array to process
                ii_s=(ii*inc_ii)-inc_ii+1;
                ii_e=(ii*inc_ii);
                jj_s=(jj*inc_jj)-inc_jj+1;
                jj_e=(jj*inc_jj);
                %Calculate indices for output array
                ind_xx=(xx*1000)-1000+ii;
                ind_yy=(yy*1000)-jj;
                
                sec_area=tempfarea(jj_s:jj_e,ii_s:ii_e);
                sec_area_mean=mean(sec_area(:));
                if ~isnan(sec_area_mean)
                    if sec_area_mean>10
                        forested_thres10(ind_xx,ind_yy)=1;
                    end
                    if sec_area_mean>50
                        forested_thres50(ind_xx,ind_yy)=1;
                    end
                    totffrac_0p01deg(ind_xx,ind_yy)=double(sum(sec_area(:))/100)/(inc_ii*inc_jj); %Total fraction covered by trees
                end
            end
            clear jj
        end
        clear ii            
        
        clear temploss tempfarea Rloss Rfarea dims_loss dims_farea
        
        fprintf('Processed %s %s. Total units is %d\n',yyind(yy,:),xxind(xx,:),cc)
    end
    clear xx
    %outtemp=['testoutthres',mat2str(yy),'.mat'];
    %save(outtemp,'forested_thres','-v7.3');
end
clear yy cc

%Calculate areas of 0.01 degree gridcells consistent with above
grid=0.01;
offset=grid/2;
lat_map=-89.995:grid:89.995;
basedist=circ_earth/(360/grid);
areag_0p01deg=NaN(length(lat_map),1);
for y=1:length(lat_map)
    areag_0p01deg(y)=basedist*basedist*cosd(lat_map(y)+offset); %m2
end
clear y basedist

%Calculate total canopy cover area
totfarea_0p01deg=totffrac_0p01deg.*repmat(areag_0p01deg',[36000 1]);

save forested_frac_0_01_degree_thres10.mat -v7.3

end %if part1


%-----
if part2

load forested_frac_0_01_degree_thres10.mat

%Now reprocess to a coarser resolution, defining a grid cell as forested
%depending on whether its forested fraction is over a threshold.
if halfdegree
    %Make calculations at 0.5 degree resolution
    nlats=360;
    nlons=720;
    lons=-179.75:0.5:179.75;
    lats=-89.75:0.5:89.75;
    conv=50;
else
    %Then do it at 1 degree resolution
    nlats=180;
    nlons=360;
    lons=-179.5:1:179.5;
    lats=-89.5:1:89.5;
    conv=100;
end

forested_thres10_regrid=zeros(nlats,nlons);
forested_thres50_regrid=zeros(nlats,nlons);
totffrac_0p01deg_regrid=zeros(nlats,nlons);
totfarea_0p01deg_regrid=zeros(nlats,nlons);
for ilath=1:nlats
    ilat_min=(ilath*conv)-conv+1;
    ilat_max=ilath*conv;
    for ilonh=1:nlons
        ilon_min=(ilonh*conv)-conv+1;
        ilon_max=ilonh*conv;
        %Find the mean open- and closed-canopy forest fraction across the grid cell
        tsec=forested_thres10(ilon_min:ilon_max,ilat_min:ilat_max);
        tsec_sum=nansum(tsec(:));
        if tsec_sum>0
            forested_thres10_regrid(ilath,ilonh)=(tsec_sum/(conv^2))*100;
        end
        clear tsec
        %Find the mean closed-canopy forest fraction across the grid cell
        tsec=forested_thres50(ilon_min:ilon_max,ilat_min:ilat_max);
        tsec_sum=nansum(tsec(:));
        if tsec_sum>0
            forested_thres50_regrid(ilath,ilonh)=(tsec_sum/(conv^2))*100;
        end
        clear tsec
        %Find the mean canopy cover fraction across the grid cell
        tsec=totffrac_0p01deg(ilon_min:ilon_max,ilat_min:ilat_max);
        tsec_sum=nansum(tsec(:));
        if tsec_sum>0
            totffrac_0p01deg_regrid(ilath,ilonh)=(tsec_sum/(conv^2))*100;
        end
        clear tsec
        %Find the total canopy area across the grid cell
        tsec=totfarea_0p01deg(ilon_min:ilon_max,ilat_min:ilat_max);
        totfarea_0p01deg_regrid(ilath,ilonh)=nansum(tsec(:));
        clear tsec
        
    end
    fprintf('ilath is %d\n',ilath)
end
clear ilath ilonh
clear ilon_min ilon_max
clear ilat_min ilat_max


%Write output to netcdf file
ncid = netcdf.create(outfile, 'NETCDF4');

dimid_lon=netcdf.defDim(ncid,'Longitude',nlons);
dimid_lat=netcdf.defDim(ncid,'Latitude',nlats);
varid_lon=netcdf.defVar(ncid,'Longitude','double',dimid_lon);
varid_lat=netcdf.defVar(ncid,'Latitude','double',dimid_lat);

netcdf.putVar(ncid,varid_lon,lons)
netcdf.putVar(ncid,varid_lat,lats)
netcdf.putAtt(ncid,varid_lon,'Units','degrees_east')
netcdf.putAtt(ncid,varid_lat,'Units','degrees_north')

varid1=netcdf.defVar(ncid,'forested_10_percent','int',[dimid_lon dimid_lat]);
netcdf.defVarDeflate(ncid,varid1,true,true,9)
netcdf.putVar(ncid,varid1,forested_thres10_regrid')
netcdf.putAtt(ncid,varid1,'Units','%')
netcdf.putAtt(ncid,varid1,'Note','Area defined as forest based on 10% canopy cover threshold in 0.01 x 0.01 degree gridcells')

varid2=netcdf.defVar(ncid,'forested_50_percent','int',[dimid_lon dimid_lat]);
netcdf.defVarDeflate(ncid,varid2,true,true,9)
netcdf.putVar(ncid,varid2,forested_thres50_regrid')
netcdf.putAtt(ncid,varid2,'Units','%')
netcdf.putAtt(ncid,varid2,'Note','Area defined as forest based on 50% canopy cover threshold in 0.01 x 0.01 degree gridcells')

varid3=netcdf.defVar(ncid,'canopy_cover_frac','float',[dimid_lon dimid_lat]);
netcdf.defVarDeflate(ncid,varid3,true,true,9)
netcdf.putVar(ncid,varid3,totffrac_0p01deg_regrid')
netcdf.putAtt(ncid,varid3,'Units','%')
netcdf.putAtt(ncid,varid3,'Note','Actual canopy cover fraction')

varid4=netcdf.defVar(ncid,'canopy_cover_area','float',[dimid_lon dimid_lat]);
netcdf.defVarDeflate(ncid,varid4,true,true,9)
netcdf.putVar(ncid,varid4,totfarea_0p01deg_regrid')
netcdf.putAtt(ncid,varid4,'Units','m2')
netcdf.putAtt(ncid,varid4,'Note','Actual canopy cover area')

glo_varid=netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,glo_varid,'Note 1','Created using Hansen et al. (2013, Science) forest cover data (year 2000) using hansen_forest_canopy_frac_calc.m')
netcdf.putAtt(ncid,glo_varid,'Note 2',date)
netcdf.putAtt(ncid,glo_varid,'Note 3','Creator: Thomas Pugh, t.a.m.pugh@bham.ac.uk')

netcdf.close(ncid)

end %if part2
