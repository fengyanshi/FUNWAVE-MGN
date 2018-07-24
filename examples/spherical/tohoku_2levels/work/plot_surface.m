clear all
fdir1='/Volumes/Seagate Backup Plus Drive/Babak/CresentCity/2arcmin/';
fdir2='/Volumes/Seagate Backup Plus Drive/Babak/CresentCity/4arcsec/';
fdir3='/Volumes/Seagate Backup Plus Drive/Babak/CresentCity/16arcsec/';

%dep=load('../external_files/depth_30min.txt');
m=4800;
n=3600;

%[n,m]=size(dep);
dx=1/30.;;
dy=1/30.;;
x=[0:m-1]*dx+132.01667;
y=[0:n-1]*dy-59.98333;


nfile=[11];
hr={'10'};

wid=8;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf


    
%fnum=sprintf('%.5d',nfile(num));
eta=load([fdir1 'eta_0011']);

%eta(eta>2)=NaN;

%subplot(length(nfile),1, num)

pcolor(x,y,eta),shading flat
hold on
caxis([-0.5 0.5])
%title([' Time = ' hr{num} ' hr '])

%caxis([-0.05 0.05])

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')
cbar=colorbar;
set(get(cbar,'ylabel'),'String','\eta (m) ')
set(gcf,'Renderer','zbuffer');

%print -djpeg eta_inlet_shoal_irr.jpg