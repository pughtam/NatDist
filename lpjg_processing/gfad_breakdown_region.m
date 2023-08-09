function [gfad_fage_reg,gfad_upper_fage_reg,gfad_lower_fage_reg]=gfad_breakdown_region(...
    gfad_file,gfad_lower_file,gfad_upper_file,...
    use_cvegmask,use_bmask,use_fmask,use_ffrac,...
    cvegmask,bmask,ffrac,fmask,rmask,nregion)
% Split up the GFAD age class data into regions.
%
% Dependencies:
% - global_grid_area.m
%
% T. A .M. Pugh
% t.a.m.pugh@bham.ac.uk
% 14.01.21

%--- Read in the GFAD data

gfad_in=ncread(gfad_file,'age');
gfad_size=size(gfad_in);
nages=gfad_size(4);
gfad_fage=squeeze(sum(gfad_in,3)); %Sum up over all 4 GFAD PFTs
clear gfad_in

gfad_upper_in=ncread(gfad_upper_file,'age');
gfad_upper_in(gfad_upper_in==-9999)=0; %Set no-data values to zero
gfad_upper_fage=squeeze(sum(gfad_upper_in,3)); %Sum up over all 4 GFAD PFTs
clear gfad_upper_in

gfad_lower_in=ncread(gfad_lower_file,'age');
gfad_lower_size=size(gfad_lower_in);
nages_lower=gfad_lower_size(4);
gfad_lower_in(gfad_lower_in==-9999)=0; %Set no-data values to zero
gfad_lower_fage=squeeze(sum(gfad_lower_in,3)); %Sum up over all 4 GFAD PFTs
clear gfad_lower_in
%Pad the lower array to the same number of age classes
gfad_lower_fage=cat(3,gfad_lower_fage(:,:,1:nages_lower-1),zeros(gfad_size(1),gfad_size(2),nages-nages_lower),gfad_lower_fage(:,:,nages_lower));
clear gfad_size gfad_lower_size nages_lower


%--- Read in forest or overall area masks ---

if use_cvegmask && use_fmask && use_bmask
    mask=cvegmask.*fmask.*bmask;
elseif use_cvegmask && use_fmask
    mask=cvegmask.*fmask;
elseif use_fmask && use_bmask
    mask=fmask.*bmask;
elseif use_cvegmask
    mask=cvegmask;
elseif use_fmask
    mask=fmask;
elseif use_bmask
    mask=bmask;
else
    mask=ones(size(cvegmask));
end
mask=fliplr(mask');
ffrac=fliplr(ffrac');

%--- Calculate regional and global age distributions ---

% Calculate unmasked grid-cell area
garea=global_grid_area()';

% Convert fractions to areas
gfad_fage_totfor=squeeze(nansum(gfad_fage,3)); %Total forest fraction
gfad_upper_fage_totfor=squeeze(nansum(gfad_fage,3)); %Total forest fraction
gfad_lower_fage_totfor=squeeze(nansum(gfad_fage,3)); %Total forest fraction
if use_ffrac
    % Standardise by total forest fraction (i.e. convert to fraction of gridcell to allow to use forest fraction from another database)
    gfad_fage_frac=gfad_fage(:,:,:)./repmat(gfad_fage_totfor(:,:),[1 1 nages]);
    gfad_upper_fage_frac=gfad_upper_fage(:,:,:)./repmat(gfad_upper_fage_totfor(:,:),[1 1 nages]);
    gfad_lower_fage_frac=gfad_lower_fage(:,:,:)./repmat(gfad_lower_fage_totfor(:,:),[1 1 nages]);
    % Convert to areas using provided forest fraction file
    gfad_fage_area=gfad_fage_frac(:,:,:).*repmat(ffrac.*mask.*garea,[1 1 nages]);
    gfad_upper_fage_area=gfad_upper_fage_frac(:,:,:).*repmat(ffrac.*mask.*garea,[1 1 nages]);
    gfad_lower_fage_area=gfad_lower_fage_frac(:,:,:).*repmat(ffrac.*mask.*garea,[1 1 nages]);
else
    gfad_fage_frac=gfad_fage;
    gfad_upper_fage_frac=gfad_upper_fage;
    gfad_lower_fage_frac=gfad_lower_fage;
    % Convert to areas directly as forest fraction is already incorporated in the data
    gfad_fage_area=gfad_fage_frac(:,:,:).*repmat(mask.*garea,[1 1 nages]);
    gfad_upper_fage_area=gfad_upper_fage_frac(:,:,:).*repmat(mask.*garea,[1 1 nages]);
    gfad_lower_fage_area=gfad_lower_fage_frac(:,:,:).*repmat(mask.*garea,[1 1 nages]);
end

% Aggregate age distributions over regions
rmask_flip=fliplr(rmask');
gfad_fage_reg=NaN(nregion,nages);
gfad_upper_fage_reg=NaN(nregion,nages);
gfad_lower_fage_reg=NaN(nregion,nages);
for rr=1:nregion
    for aa=1:nages
        gfad_fage_area_sel=squeeze(gfad_fage_area(:,:,aa));
        gfad_upper_fage_area_sel=squeeze(gfad_upper_fage_area(:,:,aa));
        gfad_lower_fage_area_sel=squeeze(gfad_lower_fage_area(:,:,aa));
        gfad_fage_reg(rr,aa)=squeeze(nansum(gfad_fage_area_sel(rmask_flip==rr)))/1e12; %Million km2
        gfad_upper_fage_reg(rr,aa)=squeeze(nansum(gfad_upper_fage_area_sel(rmask_flip==rr)))/1e12;
        gfad_lower_fage_reg(rr,aa)=squeeze(nansum(gfad_lower_fage_area_sel(rmask_flip==rr)))/1e12;
    end
end
clear aa rr gfad_fage_area_sel

