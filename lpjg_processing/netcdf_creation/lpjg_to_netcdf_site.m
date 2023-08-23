% Write netcdf files for all data underlying manuscript figures for data deposition.
% All other calculations with LPJ-GUESS data should be based on these netcdf files
% The basis is tslice files from LPJ-GUESS which give means over the years 2001-2014.
%
% Dependencies: write_netcdf_lpjg_site.m
%
% T. Pugh
% 23.08.23

addpath('./')

input_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_recovery_sims';
output_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';

sims={'boreal_fixclim_IBS150_BNE300','borealNA_fixclim_IBS150BNE300'};

dtype={'nodist','nodist'};
dreal={'site_recovery_Eurasia','site_recovery_America'};
nsims=length(sims);

PFTs={'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G'};

if (length(dtype)~=nsims || length(dreal)~=nsims)
    error('Arrays describing simulations are not consistent in length')
end

%First write files for the simulations with age forcing for the full period
for ss=1:nsims
    
    % Get C mass data    
    cmass=dlmread([input_dir,'/',sims{ss},'/cmass.out'],'',1,0);
    cmass_dim=size(cmass);

    year=cmass(:,3);
    nyear=max(year)-min(year)+1;
    nsites=length(cmass)/nyear;
    maxcmass=max(cmass(:,16));

    cmass_site=NaN(nsites,nyear,cmass_dim(2));
    for nn=1:nsites
        yy_min=nn*nyear-nyear+1;
        yy_max=nn*nyear;
        cmass_site(nn,:,:)=cmass(yy_min:yy_max,:);
    end
    clear nn yy_min yy_max

    cmass_site=cmass_site(:,:,1:length(PFTs)+3); % Get rid of the total column

    % Write to netcdf
    variable='Cveg';
    variable_longname='Carbon stock in live vegetation';
    units='kg C m-2';
    disttype=dtype{ss};
    realisation=dreal{ss};
    modname='LPJ-GUESS';
    note1='';
    
    write_netcdf_lpjg_site(cmass_site,variable,modname,disttype,realisation,...
        units,variable_longname,output_dir,note1,PFTs,'PFT');

end
clear ss
