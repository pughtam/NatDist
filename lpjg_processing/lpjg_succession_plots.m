% Make plots of biomass by PFT over a succesional sequence.
% Reads directly from LPJ-GUESS text output files, rather than netcdf
%
% T. Pugh 27.07.23

%region=2; %1 = Eurasia boreal, 2 = North America boreal

plot_years=300; % Number of years to include on the plot

lpjg_dir='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/netcdfs_for_deposition/';

for region=1:2

    if region==1
        lpjg_file='Cveg_LPJ-GUESS_site_recovery_Eurasia_nodist.nc';
        titles={'(a) Arkhangelsk region','(b) Northern Finland (Siren)','(c) Southern boreal central Russia',...
            '(d) Northern boreal central Russia'};
        pft_trans_obs1=[70,70,90,60]; % First year for transition point between broadleaved and needleleaved in observations (where only a single year is given, take +/- 10 years
        pft_trans_obs2=[90,90,110,80]; % Second year for transition point between broadleaved and needleleaved in observations
    elseif region==2
        lpjg_file='Cveg_LPJ-GUESS_site_recovery_America_nodist.nc';
        titles={'(e) Lac-Duparquet, Canada','(f) Lake Matagami Lowland, Canada','(g) Tanana River, Alaska','(h) Ontario, Canada','(i) Northern Rocky Mountains, Canada'};
        pft_trans_obs1=[140,90,100,20,140]; % First year for transition point between broadleaved and needleleaved in observations
        pft_trans_obs2=[160,110,125,40,160]; % Second year for transition point between broadleaved and needleleaved in observations
    end


    pfts={'BNE','BINE','IBS'};

    cmass=ncread([lpjg_dir,'/',lpjg_file],'Cveg'); 
    cmass_dim=size(cmass);
    nsites=cmass_dim(1);
    nyear=cmass_dim(2);
    clear cmass_dim

    maxcmass=max(cmass(:));

    if region==1
        figure
    end
    for nn=1:nsites
        if region==1
            subplot(3,3,nn)
        elseif region==2
            subplot(3,3,nn+4)
        end
        plot(1:nyear,squeeze(cmass(nn,:,[1 2 6])),'linewidth',2)
        hold on
        set(gca,'YLim',[0 maxcmass],'XLim',[0 plot_years])
        ps=polyshape([pft_trans_obs1(nn) pft_trans_obs2(nn) pft_trans_obs2(nn) pft_trans_obs1(nn)],[0 0 maxcmass maxcmass]);
        p1=plot(ps);
        set(p1,'linestyle','none','FaceColor',[0.3 0.3 0.3])
        title(titles{nn})
    end
    if region==2
        legend(pfts)
    end

end

for nn=1:9
    if nn<=6
        subplot(3,3,nn)
        set(gca,'XTickLabel','')
    else
        subplot(3,3,nn)
        xlabel('Years')
    end
    if nn==1 || nn==4 || nn==7
        subplot(3,3,nn)
        ylabel('kg C m^{-2}')
    else
        subplot(3,3,nn)
        set(gca,'YTickLabel','')
    end
end
clear nn
