% Plot age distributions from the LPJ-GUESS simulations by GFAD regions.
%
% Dependencies:
% - lpj_to_grid_func_centre
% - readmasks_func.m
% - gfad_regions.m
% - global_grid_func.m
%
% T. Pugh
% 14.01.21

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes

plotgfad=true; %Whether to also plot the GFAD data

fname='ageclass_2014';
lpjg_dir_bestest='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000/';
lpjg_dir_lowest='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/low_2se_adjparam_latosa4000/';
lpjg_dir_highest='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/high_2se_adjparam_latosa4000/';

lpjg_dir_bestest_luh='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000_luh2/';
lpjg_dir_lowest_luh='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000_luh2low/';
lpjg_dir_highest_luh='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000_luh2high/';

fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';

gfad_file='/Users/pughtam/data/GFAD_V1-1/GFAD_V1-1.nc';
gfad_upper_file='/Users/pughtam/data/GFAD_V1-1/GFAD_V1-1_upperbound.nc';
gfad_lower_file='/Users/pughtam/data/GFAD_V1-1/GFAD_V1-1_lowerbound.nc';

%---
% Read in LPJG age distributions

% Natural-only simulations
age_bestest=squeeze(lpj_to_grid_func_centre([lpjg_dir_bestest,'/',fname],1,0)); %Best estimate
s=size(age_bestest);
nage=s(3);
clear s
age_lowest=squeeze(lpj_to_grid_func_centre([lpjg_dir_lowest,'/',fname],1,0)); %Lower bound
age_highest=squeeze(lpj_to_grid_func_centre([lpjg_dir_highest,'/',fname],1,0)); %Upper bound

% LUH2 simulations
age_bestest_luh=squeeze(lpj_to_grid_func_centre([lpjg_dir_bestest_luh,'/',fname],1,0)); %Best estimate
age_lowest_luh=squeeze(lpj_to_grid_func_centre([lpjg_dir_lowest_luh,'/',fname],1,0)); %Lower bound
age_highest_luh=squeeze(lpj_to_grid_func_centre([lpjg_dir_highest_luh ,'/',fname],1,0)); %Upper bound

% Mask arrays
[cvegmask,fmask,bmask,ffrac]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir_bestest,fmask_dir,fmask_file,bmask_dir);

if use_cvegmask
    age_bestest=age_bestest.*repmat(cvegmask,[1 1 nage]);
    age_lowest=age_lowest.*repmat(cvegmask,[1 1 nage]);
    age_highest=age_highest.*repmat(cvegmask,[1 1 nage]);
    age_bestest_luh=age_bestest_luh.*repmat(cvegmask,[1 1 nage]);
    age_lowest_luh=age_lowest_luh.*repmat(cvegmask,[1 1 nage]);
    age_highest_luh=age_highest_luh.*repmat(cvegmask,[1 1 nage]);
end

if use_fmask
    age_bestest=age_bestest.*repmat(fmask,[1 1 nage]);
    age_lowest=age_lowest.*repmat(fmask,[1 1 nage]);
    age_highest=age_highest.*repmat(fmask,[1 1 nage]);
    age_bestest_luh=age_bestest_luh.*repmat(fmask,[1 1 nage]);
    age_lowest_luh=age_lowest_luh.*repmat(fmask,[1 1 nage]);
    age_highest_luh=age_highest_luh.*repmat(fmask,[1 1 nage]);
end

if use_bmask
    age_bestest=age_bestest.*repmat(bmask,[1 1 nage]);
    age_lowest=age_lowest.*repmat(bmask,[1 1 nage]);
    age_highest=age_highest.*repmat(bmask,[1 1 nage]);
    age_bestest_luh=age_bestest_luh.*repmat(bmask,[1 1 nage]);
    age_lowest_luh=age_lowest_luh.*repmat(bmask,[1 1 nage]);
    age_highest_luh=age_highest_luh.*repmat(bmask,[1 1 nage]);
end

% Aggregate to regional age class distributions based on forest cover

[rmask,regions,regions_short]=gfad_regions(gfad_file,false);
nregion=length(regions);
rmask=flipud(rmask');

garea=global_grid_area();

age_bestest_area=age_bestest.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_lowest_area=age_lowest.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_highest_area=age_highest.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_bestest_luh_area=age_bestest_luh.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_lowest_luh_area=age_lowest_luh.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_highest_luh_area=age_highest_luh.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);

age_bestest_reg=NaN(nregion,nage);
age_lowest_reg=NaN(nregion,nage);
age_highest_reg=NaN(nregion,nage);
age_bestest_luh_reg=NaN(nregion,nage);
age_lowest_luh_reg=NaN(nregion,nage);
age_highest_luh_reg=NaN(nregion,nage);
for cc=1:nage
    age_bestest_area_temp=squeeze(age_bestest_area(:,:,cc));
    age_lowest_area_temp=squeeze(age_lowest_area(:,:,cc));
    age_highest_area_temp=squeeze(age_highest_area(:,:,cc));
    age_bestest_luh_area_temp=squeeze(age_bestest_luh_area(:,:,cc));
    age_lowest_luh_area_temp=squeeze(age_lowest_luh_area(:,:,cc));
    age_highest_luh_area_temp=squeeze(age_highest_luh_area(:,:,cc));
    for nn=1:nregion
        age_bestest_reg(nn,cc)=nansum(age_bestest_area_temp(rmask==nn))/1e12; %Million km2
        age_lowest_reg(nn,cc)=nansum(age_lowest_area_temp(rmask==nn))/1e12;
        age_highest_reg(nn,cc)=nansum(age_highest_area_temp(rmask==nn))/1e12;
        age_bestest_luh_reg(nn,cc)=nansum(age_bestest_luh_area_temp(rmask==nn))/1e12;
        age_lowest_luh_reg(nn,cc)=nansum(age_lowest_luh_area_temp(rmask==nn))/1e12;
        age_highest_luh_reg(nn,cc)=nansum(age_highest_luh_area_temp(rmask==nn))/1e12;
    end
end
clear nn cc age_bestest_area_temp age_lowest_area_temp age_highest_area_temp

%---
% Optionally read the GFAD data

if plotgfad
    [gfad_fage_reg,gfad_upper_fage_reg,gfad_lower_fage_reg]=gfad_breakdown_region(...
        gfad_file,gfad_lower_file,gfad_upper_file,...
        use_cvegmask,use_bmask,use_fmask,...
        cvegmask,bmask,ffrac,rmask,nregion);
end


%---
% Bearing in mind that best estimate does not always fit within the maximum and minimum values, find upper and lower bounds
% of the overall ranges
age_bestest_reg_min=min(cat(3,age_bestest_reg,age_lowest_reg,age_highest_reg),[],3);
age_bestest_reg_max=max(cat(3,age_bestest_reg,age_lowest_reg,age_highest_reg),[],3);

age_bestest_luh_reg_min=min(cat(3,age_bestest_luh_reg,age_lowest_luh_reg,age_highest_luh_reg),[],3);
age_bestest_luh_reg_max=max(cat(3,age_bestest_luh_reg,age_lowest_luh_reg,age_highest_luh_reg),[],3);

if plotgfad
    gfad_reg_min=min(cat(3,gfad_fage_reg,gfad_lower_fage_reg,gfad_upper_fage_reg),[],3);
    gfad_reg_max=max(cat(3,gfad_fage_reg,gfad_lower_fage_reg,gfad_upper_fage_reg),[],3);
end

%---
% Make age class plots for all relevant regions

ages=5:10:150;

figure
cc=0;
%for rr=[8 10 9 7 5 6 14 13 15 2 3 4 12 11]
for rr=[10 9 7 5 6 14 13 15 2 3 4 12 11]
    cc=cc+1;
    ss(rr)=subplot(4,4,cc);
    [p1 h1 h2]=plotyy(ages(1:14),age_bestest_luh_reg(rr,1:14),ages(15),age_bestest_luh_reg(rr,15));
    hold(p1(1)); hold(p1(2));
    set(h1,'marker','.','markersize',10,'color','k')
    set(h2,'marker','.','markersize',15,'color','k','linestyle','none')
    % Other years for baseline simulation
    plot(p1(2),ages(15),age_bestest_luh_reg_max(rr,15),'k.','markersize',5,'marker','v')
    plot(p1(2),ages(15),age_bestest_luh_reg_min(rr,15),'k.','markersize',5,'marker','^')
    % Shaded area for the sensitivity studies
    poly1_y=[age_bestest_luh_reg_max(rr,1:14),fliplr(age_bestest_luh_reg_min(rr,1:14))];
    poly1_x=[ages(1:14),fliplr(ages(1:14))];
    pp1=patch(p1(1),poly1_x,poly1_y,[0.5 0.5 0.5]);
    set(pp1,'FaceAlpha',0.3,'LineStyle','none')

    h3=plot(p1(1),ages(1:14),age_bestest_reg(rr,1:14));
    h4=plot(p1(2),ages(15),age_bestest_reg(rr,15));
    plot(p1(2),ages(15),age_bestest_reg_max(rr,15),'b','markersize',5,'marker','v')
    plot(p1(2),ages(15),age_bestest_reg_min(rr,15),'b','markersize',5,'marker','^')
    set(h3,'marker','.','markersize',10,'color','b')
    set(h4,'marker','.','markersize',15,'color','b','linestyle','none')
    set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
    poly2_y=[age_bestest_reg_max(rr,1:14),fliplr(age_bestest_reg_min(rr,1:14))];
    poly2_x=[ages(1:14),fliplr(ages(1:14))];
    pp2=patch(p1(1),poly2_x,poly2_y,[0 0 0.5]);
    set(pp2,'FaceAlpha',0.3,'LineStyle','none')

    if plotgfad
        h5=plot(p1(1),ages(1:14),gfad_fage_reg(rr,1:14));
        h6=plot(p1(2),ages(15),gfad_fage_reg(rr,15));
        plot(p1(2),ages(15),gfad_reg_max(rr,15),'g','markersize',5,'marker','v')
        plot(p1(2),ages(15),gfad_reg_min(rr,15),'g','markersize',5,'marker','^')
        set(h5,'marker','.','markersize',10,'color','g')
        set(h6,'marker','.','markersize',15,'color','g','linestyle','none')
        set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
        poly3_y=[gfad_reg_max(rr,1:14),fliplr(gfad_reg_min(rr,1:14))];
        poly3_x=[ages(1:14),fliplr(ages(1:14))];
        pp3=patch(p1(1),poly3_x,poly3_y,[0 0.5 0]);
        set(pp3,'FaceAlpha',0.3,'LineStyle','none')
    end
    
    set(p1(1),'XLim',[0 155],'YLim',[0 Inf],'Box','off', 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    set(p1(2),'XLim',[0 155],'YLim',[0 Inf], 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    if cc>4
        set(p1(1),'XTick',10:10:150,'XTickLabel',{'1-10','11-20','21-30','31-40','41-50','51-60',...
            '61-70','71-80','81-90','91-100','101-110','111-120','121-130','131-140','OG'})
        set(p1(1),'XTickLabelRotation',300)
        xlabel('Age class (years)')
    else
        set(p1(1),'XTick',10:10:150,'XTickLabel','')
    end
    if cc==1 || cc==5 || cc==9 || cc==13
        ylabel(p1(1),'Young forest area (M km^{-2})')
        %ylabel(p1(1),'Young forest area (M km^{-2})')
    end
    if cc==4 || cc==8 || cc==12 || cc==16
        ylabel(p1(2),'OG forest area (M km^{-2})')
        %ylabel(p1(2),'OG forest area (M km^{-2})')
        set(get(p1(2),'Ylabel'),'Rotation',270,'VerticalAlignment','bottom')
    end
    title(regions{rr})
end
clear rr cc

%---
% Make age class plot with a simplified set of regions

[age_bestest_regsim,regions_sim,nregionsim]=reg_simplify(age_bestest_reg,nage);
[age_bestest_regsim_min]=reg_simplify(age_bestest_reg_min,nage);
[age_bestest_regsim_max]=reg_simplify(age_bestest_reg_max,nage);
[age_bestest_luh_regsim]=reg_simplify(age_bestest_luh_reg,nage);
[age_bestest_luh_regsim_min]=reg_simplify(age_bestest_luh_reg_min,nage);
[age_bestest_luh_regsim_max]=reg_simplify(age_bestest_luh_reg_max,nage);
[gfad_fage_regsim]=reg_simplify(gfad_fage_reg,nage);
[gfad_regsim_min]=reg_simplify(gfad_reg_min,nage);
[gfad_regsim_max]=reg_simplify(gfad_reg_max,nage);

figure
cc=0;
for rr=1:nregionsim
    cc=cc+1;
    ss(rr)=subplot(2,4,cc);
    [p1 h1 h2]=plotyy(ages(1:14),age_bestest_luh_regsim(rr,1:14),ages(15),age_bestest_luh_regsim(rr,15));
    hold(p1(1)); hold(p1(2));
    set(h1,'marker','.','markersize',10,'color','k')
    set(h2,'marker','.','markersize',15,'color','k','linestyle','none')
    % Other years for baseline simulation
    plot(p1(2),ages(15),age_bestest_luh_regsim_max(rr,15),'k.','markersize',5,'marker','v')
    plot(p1(2),ages(15),age_bestest_luh_regsim_min(rr,15),'k.','markersize',5,'marker','^')
    % Shaded area for the sensitivity studies
    poly1_y=[age_bestest_luh_regsim_max(rr,1:14),fliplr(age_bestest_luh_regsim_min(rr,1:14))];
    poly1_x=[ages(1:14),fliplr(ages(1:14))];
    pp1=patch(p1(1),poly1_x,poly1_y,[0.5 0.5 0.5]);
    set(pp1,'FaceAlpha',0.3,'LineStyle','none')

    h3=plot(p1(1),ages(1:14),age_bestest_regsim(rr,1:14));
    h4=plot(p1(2),ages(15),age_bestest_regsim(rr,15));
    plot(p1(2),ages(15),age_bestest_regsim_max(rr,15),'b','markersize',5,'marker','v')
    plot(p1(2),ages(15),age_bestest_regsim_min(rr,15),'b','markersize',5,'marker','^')
    set(h3,'marker','.','markersize',10,'color','b')
    set(h4,'marker','.','markersize',15,'color','b','linestyle','none')
    set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
    poly2_y=[age_bestest_regsim_max(rr,1:14),fliplr(age_bestest_regsim_min(rr,1:14))];
    poly2_x=[ages(1:14),fliplr(ages(1:14))];
    pp2=patch(p1(1),poly2_x,poly2_y,[0 0 0.5]);
    set(pp2,'FaceAlpha',0.3,'LineStyle','none')
    
    if plotgfad
        h5=plot(p1(1),ages(1:14),gfad_fage_regsim(rr,1:14));
        h6=plot(p1(2),ages(15),gfad_fage_regsim(rr,15));
        plot(p1(2),ages(15),gfad_regsim_max(rr,15),'g','markersize',5,'marker','v')
        plot(p1(2),ages(15),gfad_regsim_min(rr,15),'g','markersize',5,'marker','^')
        set(h5,'marker','.','markersize',10,'color','g')
        set(h6,'marker','.','markersize',15,'color','g','linestyle','none')
        set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
        poly3_y=[gfad_regsim_max(rr,1:14),fliplr(gfad_regsim_min(rr,1:14))];
        poly3_x=[ages(1:14),fliplr(ages(1:14))];
        pp3=patch(p1(1),poly3_x,poly3_y,[0 0.5 0]);
        set(pp3,'FaceAlpha',0.3,'LineStyle','none')
    end
    
    set(p1(1),'XLim',[0 155],'YLim',[0 Inf],'Box','off', 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    set(p1(2),'XLim',[0 155],'YLim',[0 Inf], 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    if cc>4
        set(p1(1),'XTick',10:10:150,'XTickLabel',{'1-10','11-20','21-30','31-40','41-50','51-60',...
            '61-70','71-80','81-90','91-100','101-110','111-120','121-130','131-140','OG'})
        set(p1(1),'XTickLabelRotation',300)
        xlabel('Age class (years)')
    else
        set(p1(1),'XTick',10:10:150,'XTickLabel','')
    end
    if cc==1 || cc==5 || cc==9 || cc==13
        ylabel(p1(1),'Young forest area (M km^{-2})')
        %ylabel(p1(1),'Young forest area (M km^{-2})')
    end
    if cc==4 || cc==8 || cc==12 || cc==16
        ylabel(p1(2),'OG forest area (M km^{-2})')
        %ylabel(p1(2),'OG forest area (M km^{-2})')
        set(get(p1(2),'Ylabel'),'Rotation',270,'VerticalAlignment','bottom')
    end
    title(regions_sim{rr})
end
clear rr cc


