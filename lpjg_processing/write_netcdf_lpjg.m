function write_netcdf_lpjg(data,variable,modname,disttype,realisation,year,...
            units,variable_longname,output_dir,note1,PFTs)
%Write out to netcdf file. This function to be called from lpjg_to_netcdfs.m
%
%T. Pugh
%20.10.21

fprintf('Creating output file\n')

datadim=size(data);
ndim=length(datadim);
if ndim==3
    npfts=datadim(3);
    if length(PFTs)~=npfts
        error('Number of PFTs in files does not correspond to length of PFT name array')
    end
else
    npfts=1;
end

minlat=-89.75;
maxlat=89.75;
minlon=-179.75;
maxlon=179.75;
gridspace=0.5;

latgrid=minlat:gridspace:maxlat;
longrid=minlon:gridspace:maxlon;
nlat=length(latgrid);
nlon=length(longrid);

fname=[variable,'_',modname,'_',realisation,'_',disttype,'.nc'];
outfile=[output_dir,'/',fname];
    
nccreate(outfile,'latitude','Dimensions',{'latitude',nlat})
ncwrite(outfile,'latitude',latgrid)
ncwriteatt(outfile,'latitude','units','degrees_north');
nccreate(outfile,'longitude','Dimensions',{'longitude',nlon})
ncwrite(outfile,'longitude',longrid)
ncwriteatt(outfile,'longitude','units','degrees_east');
if ndim==3
    nccreate(outfile,'PFT','Dimensions',{'PFT',npfts})
    ncwrite(outfile,'PFT',1:npfts)
    for nn=1:npfts
        ncwriteatt(outfile,'PFT',mat2str(nn),PFTs{nn});
    end
end
nccreate(outfile,'time','Dimensions',{'time',Inf})
ncwrite(outfile,'time',year)
ncwriteatt(outfile,'time','units','years');

if ndim==3
    nccreate(outfile,variable,'Dimensions',{'longitude','latitude','PFT','time'},'DeflateLevel',9)
    ncwrite(outfile,variable,permute(data,[2 1 3]))
else
    nccreate(outfile,variable,'Dimensions',{'longitude','latitude','time'},'DeflateLevel',9)
    ncwrite(outfile,variable,data')
end

ncwriteatt(outfile,variable,'longname',variable_longname)
ncwriteatt(outfile,variable,'units',units)
    
ncwriteatt(outfile,'/','DisturbanceType',disttype);
ncwriteatt(outfile,'/','Realisation',realisation);
ncwriteatt(outfile,'/','Note1',note1);
ncwriteatt(outfile,'/','Reference','This data is described in: Pugh, T.A.M., Seidl, R., Lindeskog, M., Liu, D., Chini, L.P., Senf, C., The anthropogenic imprint on temperate and boreal forest demography and carbon turnover');
ncwriteatt(outfile,'/','Institution','Lund University, Sweden');
ncwriteatt(outfile,'/','Contact','Thomas Pugh, thomas.pugh@nateko.lu.se');
ncwriteatt(outfile,'/','Version',['Version 1: ',date]);
