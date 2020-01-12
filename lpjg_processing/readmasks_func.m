function [cvegmask,fmask,bmask]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir,fmask_dir,fmask_file,bmask_dir)
% Create the masks for low vegetation C mass, low forested area or non temperate/boreal biomes
%
% T. Pugh
% 12.01.20


% Create mask to exclude grid cells where the vegetation C mass does not meet a threshold of 1 kg C m-2
if use_cvegmask
    cpool=lpj_to_grid_func_centre([lpjg_dir,'/cpool_2001_2014'],1,0);
    cveg=squeeze(cpool(:,:,1));
    cvegmask=NaN(size(cveg));
    cvegmask(cveg>1)=1;
else
    cvegmask=NaN;
end

% Create mask to exclude grid cells where at least 10% canopy is not reached
if use_fmask
    if ccmask
        fmask_var='forested_50_percent'; %Name of variable in fmask file to use
        ffrac=ncread([fmask_dir,'/',fmask_file],fmask_var)';
        fmask=NaN(size(ffrac));
        fmask(ffrac>=5)=1;
    else
        fmask_var='canopy_cover_frac'; %Name of variable in fmask file to use
        ffrac=ncread([fmask_dir,'/',fmask_file],fmask_var)';
        fmask=NaN(size(ffrac));
        fmask(ffrac>=10)=1;
    end
else
    fmask=NaN;
end

%Read in the biome mask and format to 0.5 x 0.5 degrees
if use_bmask
    bmask_temp=flipud(geotiffread([bmask_dir,'/temperate_biome_025degree.tif']));
    bmask_bor=flipud(geotiffread([bmask_dir,'/boreal_biome_025degree.tif']));
    bmask_025=zeros(size(bmask_temp));
    bmask_025(bmask_temp==1)=1;
    bmask_025(bmask_bor==1)=1;
    clear bmask_bor bmask_temp
    bmask=NaN(360,720);
    for xx=1:720
        for yy=1:360
            xx_s=(xx*2)-1;
            xx_e=xx*2;
            yy_s=(yy*2)-1;
            yy_e=yy*2;
            temp=bmask_025(yy_s:yy_e,xx_s:xx_e);
            bmask(yy,xx)=mode(temp(:));
        end
    end
    bmask(bmask==0)=NaN;
    clear xx yy xx_s xx_e yy_s yy_e    
    clear bmask_025
else
    bmask=NaN;
end
