%Script to map trait values from species level to LPJ-GUESS PFTs based on
%weighting of species by their relative prevalence across the different
%landscapes in Sommerfeld et al. (2018). It assumes that each landscape is
%of equivalent size and that the landscapes themselves, taken as a whole,
%as representative of the wider temperate forest.
%
%Units of WD (wood density) in landscape_tree_traits_complete_clean.csv are
%Units of H (species maximum height) are m
%PFT codes are TeBE=1, TeBS=2, IBS=3, TeNE/BNE=4, BINE=5, BNS=6
%
%T. Pugh
%23.08.19

%Get the mapping information for species to PFTs
fid=fopen('tree_species_traits_PFTmapping_tempbor.csv');
species_mapping=textscan(fid,'%q %q %q %q %q %d','headerlines',1,'delimiter',',');
fclose(fid);
species_all=species_mapping(1);
PFT_mapping=species_mapping{6};

%Load the landscape information
fid=fopen('species_v2.csv');
landsc_data=textscan(fid,'%d %q %d %q %d %q','headerlines',1,'delimiter',',');
fclose(fid);
landsc_spec=landsc_data(4);
landsc_perc=landsc_data{5};

%Load the trait data
fid=fopen('traits.csv');
trait_data=textscan(fid,'%q %f %f %f %d','headerlines',1,'delimiter',',');
fclose(fid);
trait_spec=trait_data(1);
trait_H=trait_data{3};
trait_WD=trait_data{4};

%Map the trait data to the landscapes
landsc_H=NaN(size(landsc_perc));
landsc_WD=NaN(size(landsc_perc));
for nn=1:length(landsc_perc)
    itrait = find(ismember(trait_spec{1}, landsc_spec{1}{nn})==1); 
    if ~isempty(itrait)
        landsc_H(nn)=trait_H(itrait);
        landsc_WD(nn)=trait_WD(itrait);
    end
end
clear nn itrait


nspec=length(species_mapping{1});
npft=max(PFT_mapping);

%Sum prevalence values for all species across all sites
spec_preval_sum=NaN(nspec,1);
spec_WD=NaN(nspec,1);
spec_WDstd=NaN(nspec,1);
spec_H=NaN(nspec,1);
spec_Hstd=NaN(nspec,1);
for nn=1:nspec
    ispec = find(ismember(landsc_spec{1}, species_all{1}{nn})==1);    
    spec_preval_sum(nn)=sum(landsc_perc(ispec));
    %Also extract the species-level values of the traits
    spec_WD(nn)=mean(landsc_WD(ispec));
    spec_WDstd(nn)=std(landsc_WD(ispec));
    spec_H(nn)=mean(landsc_H(ispec));
    spec_Hstd(nn)=std(landsc_H(ispec));
end
clear nn ispec

%Check that there is no intra-species variation in traits (should not be)
if spec_WDstd~=0
    error('Intra-species variation in WD')
end
if spec_Hstd~=0
    error('Intra-species variation in H')
end
clear spec_WDstd spec_Hstd

%For each PFT calculate averages of the traits across all species,
%weighted by the prevalence values
WD_PFT=NaN(npft,1);
H_PFT=NaN(npft,1);
WD_PFT_std=NaN(npft,1);
H_PFT_std=NaN(npft,1);
for pp=1:npft
    ipft=find(PFT_mapping==pp & isnan(spec_WD)==0 & isnan(spec_H)==0);
    WD_PFT(pp)=wmean(spec_WD(ipft),spec_preval_sum(ipft));
    H_PFT(pp)=wmean(spec_H(ipft),spec_preval_sum(ipft));
    WD_PFT_std(pp)=sqrt(var(spec_WD(ipft),spec_preval_sum(ipft)));
    H_PFT_std(pp)=sqrt(var(spec_H(ipft),spec_preval_sum(ipft)));
end
clear pp

%Convert wood density units from mg mm-3 to kgC m-3
DM_to_C=0.5;
mg_to_kg=1e-6;
mm3_to_m3=1e9;
WD_PFT=WD_PFT*DM_to_C*mg_to_kg*mm3_to_m3;
WD_PFT_std=WD_PFT_std*DM_to_C*mg_to_kg*mm3_to_m3;

pfts={'TeBE','TeBS','IBS','TeNE/BNE','BINE','BNS'};
figure
subplot(2,1,1)
hold on
bar(WD_PFT)
e1=errorbar(1:6,WD_PFT,WD_PFT_std);
set(e1,'linestyle','none')
set(gca,'XTick',1:6,'XTickLabel',pfts)
ylabel('Wood density (kg C m^{-3})')

subplot(2,1,2)
hold on
bar(H_PFT)
e2=errorbar(1:6,H_PFT,H_PFT_std);
set(e2,'linestyle','none')
set(gca,'XTickLabel',pfts)
set(gca,'XTick',1:6,'XTickLabel',pfts)
ylabel('Maximum height (m)')

fprintf('Wood density\n')
fprintf('%s %s %s %s %s %s\n','TeBE','TeBS','IBS','TeNE/BNE','BINE','BNS')
fprintf('%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',WD_PFT);
fprintf('%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',WD_PFT_std);

fprintf('Maximum height\n')
fprintf('%s %s %s %s %s %s\n','TeBE','TeBS','IBS','TeNE/BNE','BINE','BNS')
fprintf('%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',H_PFT);
fprintf('%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',H_PFT_std);
