

use_cvegmask=true; %Mask by a minimum simulated vegetation biomass density
use_fmask=true; %Mask by current forest area
ccmask=false; %Use a closed-canopy forest mask (if use_fmask=true)
use_bmask=true; %Mask by temperate/boreal biomes


fname='ageclass_2014';
lpjg_dir_bestest='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/best_est_adjparam_latosa4000/';
lpjg_dir_lowest='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/low_2se_adjparam_latosa4000/';
lpjg_dir_highest='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/LPJG_results/high_2se_adjparam_latosa4000/';


fmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal';
fmask_file='hansen_forested_canopy_frac_0p5deg.nc4';
bmask_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/biomes/From_Cornelius_inc_boreal';
regmaskfile_gfad='/Users/pughtam/data/GFAD_V1-1/GFAD_V1-1.nc';



% Read in LPJG age distributions

% Best estimate
age_bestest=squeeze(lpj_to_grid_func_centre([lpjg_dir_bestest,'/',fname],1,0));
s=size(age_bestest);
nage=s(3);
clear s

% Lower bound
age_lowest=squeeze(lpj_to_grid_func_centre([lpjg_dir_lowest,'/',fname],1,0));

% Upper bound
age_highest=squeeze(lpj_to_grid_func_centre([lpjg_dir_highest,'/',fname],1,0));


% Mask arrays
[cvegmask,fmask,bmask,ffrac]=readmasks_func(use_cvegmask,use_fmask,ccmask,use_bmask,lpjg_dir_bestest,fmask_dir,fmask_file,bmask_dir);

if use_cvegmask
    age_bestest=age_bestest.*repmat(cvegmask,[1 1 nage]);
    age_lowest=age_lowest.*repmat(cvegmask,[1 1 nage]);
    age_highest=age_highest.*repmat(cvegmask,[1 1 nage]);
end

if use_fmask
    age_bestest=age_bestest.*repmat(fmask,[1 1 nage]);
    age_lowest=age_lowest.*repmat(fmask,[1 1 nage]);
    age_highest=age_highest.*repmat(fmask,[1 1 nage]);
end

if use_bmask
    age_bestest=age_bestest.*repmat(bmask,[1 1 nage]);
    age_lowest=age_lowest.*repmat(bmask,[1 1 nage]);
    age_highest=age_highest.*repmat(bmask,[1 1 nage]);
end

% Aggregate to regional age class distributions based on forest cover

onedeg=false;
[rmask,regions,regions_short]=gfad_regions(regmaskfile_gfad,onedeg);
clear onedeg
nregion=length(regions);
rmask=flipud(rmask');

garea=global_grid_area();

age_bestest_area=age_bestest.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_lowest_area=age_lowest.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);
age_highest_area=age_highest.*repmat(ffrac,[1,1,nage]).*repmat(garea,[1,1,nage]);

age_bestest_reg=NaN(nregion,nage);
age_lowest_reg=NaN(nregion,nage);
age_highest_reg=NaN(nregion,nage);
for cc=1:nage
    age_bestest_area_temp=squeeze(age_bestest_area(:,:,cc));
    age_lowest_area_temp=squeeze(age_lowest_area(:,:,cc));
    age_highest_area_temp=squeeze(age_highest_area(:,:,cc));
    for nn=1:nregion
        age_bestest_reg(nn,cc)=nansum(age_bestest_area_temp(rmask==nn));
        age_lowest_reg(nn,cc)=nansum(age_lowest_area_temp(rmask==nn));
        age_highest_reg(nn,cc)=nansum(age_highest_area_temp(rmask==nn));
    end
end
clear nn cc age_bestest_area_temp age_lowest_area_temp age_highest_area_temp

% Bearing in mind that best estimate does not always fit within the maximum and minimum values, find upper and lower bounds
% of the overall ranges
age_bestest_reg_min=min(cat(3,age_bestest_reg,age_lowest_reg,age_highest_reg),[],3);
age_bestest_reg_max=max(cat(3,age_bestest_reg,age_lowest_reg,age_highest_reg),[],3);

% Make age class plots


ages=5:10:150;

%ind1=find(dist_year==2015);

figure
cc=0;
%for rr=[8 10 9 7 5 6 14 13 15 2 3 4 12 11]
for rr=[10 9 7 5 6 14 13 15]
    cc=cc+1;
    ss(rr)=subplot(2,4,cc);
    [p1 h1 h2]=plotyy(ages(1:14),age_bestrest_reg(rr,1:14),ages(15),age_bestrest_reg(rr,15));
    hold(p1(1)); hold(p1(2));
    set(h1,'marker','.','markersize',10,'color','k')
    set(h2,'marker','.','markersize',15,'color','k','linestyle','none')
    % Other years for baseline simulation
    plot(p1(2),ages(15),age_bestest_reg_max(rr,15),'k.','markersize',5,'marker','v')
    plot(p1(2),ages(15),age_bestest_reg_min(rr,15),'k.','markersize',5,'marker','^')
    % Shaded area for the sensitivity studies
    poly1_y=[age_bestest_reg_max(rr,1:14),fliplr(age_bestest_reg_min(rr,1:14))];
    poly1_x=[ages(1:14),fliplr(ages(1:14))];
    pp1=patch(p1(1),poly1_x,poly1_y,[0.5 0.5 0.5]);
    set(pp1,'FaceAlpha',0.3,'LineStyle','none')

%     if plotextraline
%         h3=plot(p1(1),ages(1:14),distplot2_baseline(ind1,1:14,rr));
%         h4=plot(p1(2),ages(15),distplot2_baseline(ind1,15,rr));
%         plot(p1(2),ages(15),distplot2_baseline_upperrange(ind1,15,rr),'b','markersize',5,'marker','v')
%         plot(p1(2),ages(15),distplot2_baseline_lowerrange(ind1,15,rr),'b','markersize',5,'marker','^')
%         set(h3,'marker','.','markersize',10,'color','b')
%         set(h4,'marker','.','markersize',15,'color','b','linestyle','none')
%         set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
%         poly2_y=[distplot2_baseline_upperrange(ind1,1:14,rr),fliplr(distplot2_baseline_lowerrange(ind1,1:14,rr))];
%         poly2_x=[ages(1:14),fliplr(ages(1:14))];
%         pp2=patch(p1(1),poly2_x,poly2_y,[0 0 0.5]);
%         set(pp2,'FaceAlpha',0.3,'LineStyle','none')
%     end
% 
%     h5=plot(p1(1),ages(1:14),gfadplot(1:14,rr));
%     h6=plot(p1(2),ages(15),gfadplot(15,rr));
%     plot(p1(2),ages(15),gfadplot_upperrange(15,rr),'g','markersize',5,'marker','v')
%     plot(p1(2),ages(15),gfadplot_lowerrange(15,rr),'g','markersize',5,'marker','^')
%     set(h5,'marker','.','markersize',10,'color','g')
%     set(h6,'marker','.','markersize',15,'color','g','linestyle','none')
%     set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
%     poly3_y=[gfadplot_upperrange(1:14,rr)',fliplr(gfadplot_lowerrange(1:14,rr)')];
%     poly3_x=[ages(1:14),fliplr(ages(1:14))];
%     pp3=patch(p1(1),poly3_x,poly3_y,[0 0.5 0]);
%     set(pp3,'FaceAlpha',0.3,'LineStyle','none')
    
    
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
        ylabel(p1(1),'CC young forest area (M km^{-2})')
        %ylabel(p1(1),'Young forest area (M km^{-2})')
    end
    if cc==4 || cc==8 || cc==12 || cc==16
        ylabel(p1(2),'CC OG forest area (M km^{-2})')
        %ylabel(p1(2),'OG forest area (M km^{-2})')
        set(get(p1(2),'Ylabel'),'Rotation',270,'VerticalAlignment','bottom')
    end
    title(regions{rr})
end
clear rr cc


