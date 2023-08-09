% Extract disturbance retrun periods from LPJ-GUESS simulations for specific sites to compare with the literature
%
% T. Pugh
% 09.08.23

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';
lpjg_file='distprob_LPJ-GUESS_standard_nat_2014.nc';

dist=permute(ncread([lpjg_dir,'/',lpjg_file],'distprob'),[2 1]); 
dist(dist==0)=NaN;
lat=ncread([lpjg_dir,'/',lpjg_file],'latitude');
lon=ncread([lpjg_dir,'/',lpjg_file],'longitude');

%% Now extract the disturbance return periods for sites in the literature

% Kelly et al. (2013), 111 year return period
indlat=find(lat==65.75);
indlon=find(lon==-146.25);

fprintf('Kelly et al. (2013): %f\n',1./dist(indlat,indlon))

% Kharuk et al. (2016), site I, 106 +/- 36
indlat=find(lat==65.25);
indlon=find(lon==99.75);

fprintf('Kharuk et al. (2016) Site I: %f\n',1./dist(indlat,indlon))

% Kharuk et al. (2016), site II, 82 +/- 7
indlat=find(lat==63.75);
indlon=find(lon==105.75);

fprintf('Kharuk et al. (2016) Site II: %f\n',1./dist(indlat,indlon))

% Kharuk et al. (2016), site III, 200 +/- 51
indlat=find(lat==66.25);
indlon=find(lon==99.75);

fprintf('Kharuk et al. (2016) Site III: %f\n',1./dist(indlat,indlon))

% Prince et al. (2018), 120 (100-142, 95% confidence interval)
indlat=find(lat==60.25);
indlon=find(lon==-134.75);

fprintf('Prince et al. (2018): %f\n',1./dist(indlat,indlon))

% Senichi et al. (2015), 200
indlat=find(lat==49.25);
indlon=find(lon==-89.75);

fprintf('Senichi et al. (2015): %f\n',1./dist(indlat,indlon))

% Higuera et al. (2009), Code/Wild Tussock, 135 (113-160) year return period
indlat=find(lat==67.25);
indlon=find(lon==-151.25);

fprintf('Higuera et al. (2009): %f\n',1./dist(indlat,indlon))

% Higuera et al. (2009), Ruppert, 171 (135-216) year return period
indlat=find(lat==66.75);
indlon=find(lon==-154.25);

fprintf('Higuera et al. (2009): %f\n',1./dist(indlat,indlon))


