function [data_out,regionnames,nregionsim]=reg_simplify(data_in,nage)

nregionsim=8;

regionnames={'Canada','West US','North-East US','South-East US','Europe','West Russia','Central/East Russia','East Asia'};

data_out=NaN(nregionsim,nage);
data_out(1,:)=data_in(10,:)+data_in(9,:); %Canada
data_out(2,:)=data_in(7,:); %West US
data_out(3,:)=data_in(5,:); %North-East US
data_out(4,:)=data_in(6,:); %South-East US
data_out(5,:)=data_in(14,:)+data_in(13,:); %Europe
data_out(6,:)=data_in(2,:); %West Russia
data_out(7,:)=data_in(3,:)+data_in(4,:); %Central/East Russia
data_out(8,:)=data_in(11,:)+data_in(12,:)+data_in(15,:); %China, Japan and Korea