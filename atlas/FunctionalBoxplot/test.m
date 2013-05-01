close all;
clear all;

%default values

factor=1.5; 
fullout='False'; 
barcol='b'; 
outliercol='r'; 
color=[[0.8784, 0.6902, 1.0]; [1, 0, 1]; [0.5451, 0, 0.8] ]; 
prob=[0.75, 0.5, 0.25]; 
show='True';
method='MBD'; 
depth=[]; 

outputPrefix = './Pic';

ncasem = 39;
ncasef = 54;
nage   = 31;

fid = fopen('hgtm.dat','rt');
hgtmmat = reshape(fscanf(fid,'%f'),[nage,ncasem]);

fid = fopen('hgtf.dat','rt');
hgtfmat = reshape(fscanf(fid,'%f'),[nage,ncasef]);

age = [ 1:0.25:2, 3:8, 8.5:0.5:18 ]';

%fbplot of boys' height
figure, fbplot(hgtmmat,age,depth,method,show,prob,color,outliercol,barcol,fullout,factor)
xlim([0.5,18.5])
ylim([60,200])
xlabel('Age (Years)')
ylabel('Height (cm)')
title('Boys')
filename = sprintf( '%s/Boys.png', outputPrefix ); 
saveas( gca, filename );

%fbplot of girls' height
figure, fbplot(hgtfmat,age,depth,method,show,prob,color,outliercol,barcol,fullout,factor)
xlim([0.5,18.5])
ylim([60,200])
xlabel('Age (Years)')
ylabel('Height (cm)')
title('Girls')
filename = sprintf( '%s/Girls.png', outputPrefix );
saveas( gca, filename );
