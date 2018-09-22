clear all;  clc;  close all;

dep=load('U_ny3600_mx4800.txt');
dx=2/60;   %0.0333;
dy=2/60;   %0.0333;

[n, m] = size(dep);

skip=2;
dep4=dep(1:skip:n,1:skip:m);

clear n m;
[n, m] = size(dep4);
fname=[ 'U_4min_ny',num2str(n),'mx',num2str(m),'.dat' ];
[FileID]=fopen( fname, 'wt' );
for line=1:1:size(dep4,1);
    fprintf( FileID, '%12.4f', dep4(line,:) );
    fprintf( FileID,'\n');
end; 
clear FileID
fclose all;

x=[0:m-1]*(dx*2)+132.01667;
y=[0:n-1]*(dy*2)-59.98333;

figure(1)
pcolor(x,y,dep4),shading interp
colormap jet
colorbar
axis equal;   axis tight;
xlabel(' Lon (deg) ');   ylabel(' Lat (deg) ');

pname = ['surf_U.png'];
set(gcf,'PaperPositionMode', 'auto');
print ('-dpng', pname);