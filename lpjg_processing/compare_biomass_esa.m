% Compare LPJ-GUESS biomass with ESA CCI biomass for each of the landscapes
% 
% T. Pugh 23.07.23

calc_indata=false; % true = calculate which cells are in the landscape (first time, slow), false = read this information from precreated mat files
lpjg_from_netcdf=true; % Read LPJ-GUESS data from netcdf files or directly from testruns for the individual landscapes

% Set the data directories
esa_dir='/Users/pughtam/data/ESA_CCI_Biomass';
land_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/landscape_outlines';
lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_tests/';

% Get the ESA CCI biomass data
% Downloaded from https://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v3.0/netcdf
lon=ncread([esa_dir,'/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.nc'],'lon');
lat=ncread([esa_dir,'/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.nc'],'lat');
% Size of grid cells
%loninc=lon(2)-lon(1);
%latinc=lat(1)-lat(2);

% Get the LPJ-GUESS biomass data
if lpjg_from_netcdf
    lpjg_veg=permute(ncread([lpjg_dir,'/Cveg_LPJ-GUESS_standardIBS150BNE300normalkl_nat_2014.nc'],'Cveg'),[2 1]);
    lpjg_lon=permute(ncread([lpjg_dir,'/Cveg_LPJ-GUESS_standardIBS150BNE300normalkl_nat_2014.nc'],'longitude'),[2 1]);
    lpjg_lat=permute(ncread([lpjg_dir,'/Cveg_LPJ-GUESS_standardIBS150BNE300normalkl_nat_2014.nc'],'latitude'),[2 1]);
else
    lpjg_in=readtable('/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_site_sims/landscapes_IBS150_BNE300/cpool.out','FileType','delimitedtext');
    [~,ia,ic]=unique(lpjg_in(:,1:2));
    ia_unsorted=sort(ia); % Reverse the sorting of the indexes to preserve the original site order form the simulation
    ic_unsorted=sort(ic);
    lpjg_lon=(floor(lpjg_in.Lon(ia_unsorted)*2)/2)+0.25;
    lpjg_lat=(floor(lpjg_in.Lat(ia_unsorted)*2)/2)+0.25;
    clear lpjg_lonlat
    lpjg_veg=NaN(max(ic),1);
    for nn=1:max(ic)
        temp=lpjg_in(ic_unsorted==nn,:);
        lpjg_veg(nn)=mean(temp.VegC(temp.Year>2000 & temp.Year<=2014));
    end
    clear lpjg_in ia ic ia_unsorted ic_unsorted nn
end

[lpjg_lon_grid,lpjg_lat_grid]=meshgrid(lpjg_lon,lpjg_lat);

% Get the landscape outlines
Sb=shaperead([land_dir,'/sites_boreal.shp']);
St=shaperead([land_dir,'/sites_temperate.shp']);
S=[Sb;St];
clear Sb St
nS=length(S);

% Now iterate through the landscapes and calculate the biomass in the ESA dataset and the LPJ-GUESS simulations
esa_biomass_mean=NaN(nS,1);
esa_biomass_nanmean=NaN(nS,1);
esa_biomass_median=NaN(nS,1);
esa_biomass_sd=NaN(nS,1);
lpjg_biomass_mean=NaN(nS,1);
lpjg_biomass_std=NaN(nS,1);
for ss=1:nS

    % Get the range of the landscape in lat and lon to extract the appropriate section from the ESA biomass dataset
    xmin=min(S(ss).X);
    xmax=max(S(ss).X);
    ymin=min(S(ss).Y);
    ymax=max(S(ss).Y);

    ilonmin=find(abs(lon-xmin)==min(abs(lon-xmin)));
    ilonmax=find(abs(lon-xmax)==min(abs(lon-xmax)));
    ilatmin=find(abs(lat-ymin)==min(abs(lat-ymin)));
    ilatmax=find(abs(lat-ymax)==min(abs(lat-ymax)));
    clat=ilatmin-ilatmax;
    clon=ilonmax-ilonmin;
    nstart=[ilonmin ilatmax 1]; % Note that the latitude variable is flipped
    ncount=[clon clat 1];

    biomass_ext=ncread([esa_dir,'/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.nc'],'agb',nstart,ncount);
    lon_ext=ncread([esa_dir,'/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.nc'],'lon',ilonmin,clon);
    lat_ext=ncread([esa_dir,'/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.nc'],'lat',ilatmax,clat);
    [lon_ext_grid,lat_ext_grid]=meshgrid(lon_ext,lat_ext);

    % Find the parts of the extracted ESA biomass data that are in the landscape polygon and calculate stats
    if calc_indata
        in_esa = inpolygon(lon_ext_grid,lat_ext_grid,S(ss).X,S(ss).Y);
        save(['in_esa_landscape_',mat2str(ss),'.mat'],'in_esa','ss','nstart','ncount','lon_ext','lat_ext')
    else
        load(['in_esa_landscape_',mat2str(ss),'.mat'])
    end

    biomass_ext_nan=biomass_ext;
    biomass_ext_nan(biomass_ext==0)=NaN;

    esa_biomass_mean(ss)=mean(biomass_ext(in_esa));
    esa_biomass_nanmean(ss)=nanmean(biomass_ext_nan(in_esa));
    esa_biomass_median(ss)=median(biomass_ext(in_esa));
    esa_biomass_sd(ss)=std(biomass_ext(in_esa));

    % Now identify the appropriate grid cells in the LPJ-GUESS simulation
    if lpjg_from_netcdf
        lonr=(floor(S(ss).X*2)/2)+0.25;
        latr=(floor(S(ss).Y*2)/2)+0.25;

        ir=unique(table(latr',lonr'));

        ilpjg=NaN(height(ir),1);
        for nn=1:height(ir)
            if ~isnan(ir.Var1(nn)) && ~isnan(ir.Var2(nn))
                ilpjg(nn)=find(lpjg_lon_grid==ir.Var2(nn) & lpjg_lat_grid==ir.Var1(nn));
            end
        end
        clear nn
        ilpjg(isnan(ilpjg))=[];

        lpjg_biomass_mean(ss)=mean(lpjg_veg(ilpjg));
        lpjg_biomass_std(ss)=std(lpjg_veg(ilpjg)); % Note that this standard deviation is not comparable to that of the ESA biomass data because of the different scales that it is calculated over
    else
        % Assume that the lpj-guess simulations are conducted in the same order

        if ss==1
            % Make a diagnostic plot of lat and lon to confirm this
            landlons=NaN(length(S),1);
            landlats=NaN(length(S),1);
            for nn=1:nS
                landlons(nn)=nanmean(S(nn).X);
                landlats(nn)=nanmean(S(nn).Y);
            end
            clear nn
            figure
            subplot(1,2,1)
            plot(landlons,lpjg_lon,'.')
            title('Longitudes')
            subplot(1,2,2)
            plot(landlats,lpjg_lat,'.')
            title('Latitudes')
            fprintf('NOTE: Check that longitudes and latitudes correspond in the diagnostic plot\n')
        end

        % Then simply assign the results from the simulations
        lpjg_biomass_mean(ss)=lpjg_veg(ss);
        lpjg_biomass_std(ss)=0;
    end


    % Now convert the units for the LPJ-GUESS output to be consistent with the ESA values
    % LPJ-GUESS is in kg C m-2 total biomass
    % ESA is in Mg DM ha-1 AGB (This is defined as the mass, expressed as oven-dry weight of the woody parts (stem,
    % bark, branches and twigs) of all living trees excluding stump and roots)
    m2_per_ha=1e4;
    kg_per_Mg=1000;
    C_per_DM=0.5;
    agb_frac=0.75; %NOTE: Need to get a reference for an appropriate value

    lpjg_biomass_mean(ss)=lpjg_biomass_mean(ss)*m2_per_ha/kg_per_Mg/C_per_DM*agb_frac;
    lpjg_biomass_std(ss)=lpjg_biomass_std(ss)*m2_per_ha/kg_per_Mg/C_per_DM*agb_frac;

    fprintf('Completed %d out of %d landscapes\n',ss,nS)

end
clear ss

save('biomass_esa_lpjg.mat','lpjg_biomass_mean','lpjg_biomass_std','esa_biomass_mean','esa_biomass_sd')

figure
plot(esa_biomass_nanmean(1:49),lpjg_biomass_mean(1:49),'ro')
hold on
plot(esa_biomass_nanmean(50:77),lpjg_biomass_mean(50:77),'bo')
set(gca,'XLim',[0 350],'YLim',[0 350])
l1=line([0 350],[0 350]);
set(l1,'linestyle',':','color','k')
xlabel('ESA AGB (Mg ha^{-1})')
ylabel('LPJ-GUESS AGB (Mg ha^{-1})')


% Make a gridlist for test these sites with LPJ-GUESS
landlons=NaN(length(S),1);
landlats=NaN(length(S),1);
for nn=1:nS
    landlons(nn)=nanmean(S(nn).X);
    landlats(nn)=nanmean(S(nn).Y);
end
clear nn
writetable(table(landlons,landlats),'gridlist_landscapes.txt','delimiter',' ')


% Some checking for where overestimates are occurring
aa=find(lpjg_biomass_mean>esa_biomass_mean);
landlons=NaN(length(S),1);
landlats=NaN(length(S),1);
for nn=1:nS
    landlons(nn)=nanmean(S(nn).X);
    landlats(nn)=nanmean(S(nn).Y);
end
clear nn
figure
plot(landlons,landlats,'ko')
hold on
plot(landlons(aa),landlats(aa),'ro')
load coastlines
plot(coastlon, coastlat, 'k');


