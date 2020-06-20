% Create an input file for LPJ-GUESS that specifies whether the forest type in each grid cell is broadleaf, needleleaf or
% mixed. Assignment based on ESA landcover.
%
% Dependencies:
% - esa_broadleaf_frac_0p5deg.mat (from esa_broadleaf_frac_0p5deg.m)
%
% T. Pugh
% 20.06.20

load esa_broadleaf_frac_0p5deg.mat

lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;

% Assign grid cells to forest types:
% >80% broadleaf = broadleaf
% <20% broadleaf = needleleaf
% 20-80% broadleaf = mixed
broad=zeros(size(broadleaf_frac_esa));
broad(broadleaf_frac_esa>80)=1;
broad(isnan(broadleaf_frac_esa))=NaN;

needle=zeros(size(broadleaf_frac_esa));
needle(broadleaf_frac_esa<20)=1;
needle(isnan(broadleaf_frac_esa))=NaN;

mixed=zeros(size(broadleaf_frac_esa));
mixed(broadleaf_frac_esa>=20 & broadleaf_frac_esa<=80)=1;
mixed(isnan(broadleaf_frac_esa))=NaN;

% Format into array for output
outarray=NaN(1,5);
cc=0;
for xx=1:720
    for yy=1:360
        if isnan(broadleaf_frac_esa(yy,xx))==0
            cc=cc+1;
            outarray(cc,1)=lons(xx);
            outarray(cc,2)=lats(yy);
            outarray(cc,3)=broad(yy,xx);
            outarray(cc,4)=needle(yy,xx);
            outarray(cc,5)=mixed(yy,xx);
        end
    end
end
clear xx yy

outtable=array2table(outarray,'VariableNames',{'Lon','Lat','Broadl','Needle','Mixed'});
        
% Write to text file
writetable(outtable,'lpjg_foresttype_input.txt','Delimiter',' ')
