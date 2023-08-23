function write_netcdf_lpjg_site(data,variable,modname,disttype,realisation,...
            units,variable_longname,output_dir,note1,dim3labels,dim3name)
% Write out to netcdf file. This function to be called from lpjg_to_netcdfs_site.m
%
% T. Pugh
% 23.08.23

fprintf('Creating output file\n')

datadim=size(data);

ndim3=datadim(3);
if length(dim3labels)+3~=ndim3
    error('Number of PFTs in files does not correspond to length of PFT name array')
end


fname=[variable,'_',modname,'_',realisation,'_',disttype,'.nc'];
outfile=[output_dir,'/',fname];
    
nccreate(outfile,'latitude','Dimensions',{'nsite',datadim(1)})
ncwrite(outfile,'latitude',data(:,1,2))
ncwriteatt(outfile,'latitude','units','degrees_north');
nccreate(outfile,'longitude','Dimensions',{'nsite',datadim(1)})
ncwrite(outfile,'longitude',data(:,1,1))
ncwriteatt(outfile,'longitude','units','degrees_east');

nccreate(outfile,dim3name,'Dimensions',{dim3name,ndim3-3})
ncwrite(outfile,dim3name,1:ndim3-3)
for nn=1:ndim3-3
    ncwriteatt(outfile,dim3name,mat2str(nn),dim3labels{nn});
end

nccreate(outfile,'time','Dimensions',{'time',datadim(2)})
ncwrite(outfile,'time',data(1,:,3))
ncwriteatt(outfile,'time','units','years');

nccreate(outfile,variable,'Dimensions',{'nsite','time',dim3name},'DeflateLevel',9)
ncwrite(outfile,variable,data(:,:,4:ndim3))
ncwriteatt(outfile,variable,'longname',variable_longname)
ncwriteatt(outfile,variable,'units',units)
    
ncwriteatt(outfile,'/','DisturbanceType',disttype);
ncwriteatt(outfile,'/','Realisation',realisation);
ncwriteatt(outfile,'/','Note1',note1);
ncwriteatt(outfile,'/','Reference','This data is described in: Pugh, T.A.M., Seidl, R., Lindeskog, M., Liu, D., Chini, L.P., Senf, C., The anthropogenic imprint on temperate and boreal forest demography and carbon turnover');
ncwriteatt(outfile,'/','Institution','Lund University, Sweden');
ncwriteatt(outfile,'/','Contact','Thomas Pugh, thomas.pugh@nateko.lu.se');
ncwriteatt(outfile,'/','Version',['Version 1: ',date]);
