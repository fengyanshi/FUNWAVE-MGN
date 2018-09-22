clear all
fdir='../result/';

Nghost=4;
totalres=10000;

dep1=load('../dep_4min_ny1800mx2400.dat');
dep2=load([fdir 'Grd02_dep.out']);
dep3=load([fdir 'Grd03_dep.out']);
dep4=load([fdir 'Grd04_dep.out']);
dep5=load([fdir 'Grd05_dep.out']);

grid=load(['./subgrid_info.txt']);

[n1,m1]=size(dep1);
dx1=4/60;
dy1=4/60;
x1=[0:m1-1]*dx1+132.01667;
y1=[0:n1-1]*dy1-59.98333;
[Xgrid1, Ygrid1]=meshgrid(x1,y1);
%-------------------------------------------------------------

m2=grid(1,1);
n2=grid(1,2);
isk2=grid(1,3);
trackx2=grid(1,4);
tracky2=grid(1,5);

dx2=dx1/isk2;
dy2=dy1/isk2;

x2 = 0 : dx2 : (m2-1)*dx2;
y2 = 0 : dy2 : (n2-1)*dy2;
x2 = x2+132.01667;
y2 = y2-59.98333;

x2 = x2 + (trackx2-Nghost-1)*dx1 + (Nghost*dx2);
y2 = y2 + (tracky2-Nghost-1)*dy1 + (Nghost*dy2);
[Xgrid2, Ygrid2]=meshgrid(x2,y2);
%-------------------------------------------------------------

m3=grid(2,1);
n3=grid(2,2);
isk3=grid(2,3);
trackx3=grid(2,4);
tracky3=grid(2,5);
% 
dx3=dx2/isk3;
dy3=dy2/isk3;
% 
x3 = 0 : dx3 : (m3-1)*dx3;
y3 = 0 : dy3 : (n3-1)*dy3;
% 
% x3 = x2( trackx3(1) ) + x3;
% y3 = y2( tracky3(1) ) + y3;
x3 = x3 + min( x2 ) + (trackx3-Nghost-1)*dx2 + (Nghost*dx3);
y3 = y3 + min( y2 ) + (tracky3-Nghost-1)*dy2 + (Nghost*dy3);
[Xgrid3, Ygrid3]=meshgrid(x3,y3);
%-------------------------------------------------------------

m4=grid(3,1);
n4=grid(3,2);
isk4=grid(3,3);
trackx4=grid(3,4);
tracky4=grid(3,5);
% 
dx4=dx3/isk4;
dy4=dy3/isk4;
% 
x4 = 0 : dx4 : (m4-1)*dx4;
y4 = 0 : dy4 : (n4-1)*dy4;
% 
x4 = x4 + min( x3 ) + (trackx4-Nghost-1)*dx3 + (Nghost*dx4);
y4 = y4 + min( y3 ) + (tracky4-Nghost-1)*dy3 + (Nghost*dy4);
[Xgrid4, Ygrid4]=meshgrid(x4,y4);
%-------------------------------------------------------------

m5=grid(4,1);
n5=grid(4,2);
isk5=grid(4,3);
trackx5=grid(4,4);
tracky5=grid(4,5);
% 
dx5=dx4/isk5;
dy5=dy4/isk5;
% 
x5 = 0 : dx5 : (m5-1)*dx5;
y5 = 0 : dy5 : (n5-1)*dy5;
% 
x5 = x5 + min( x4 ) + (trackx5-Nghost-1)*dx4 + (Nghost*dx5);
y5 = y5 + min( y4 ) + (tracky5-Nghost-1)*dy4 + (Nghost*dy5);
[Xgrid5, Ygrid5]=meshgrid(x5,y5);
%-------------------------------------------------------------

% nfile=[2 9];
% hr={'1' '8'};

% wid=5;
% len=7;
% set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
% clf

time=load([fdir 'Grd01_track.txt']);
time=time/3600;

for num=0:totalres;
    
fnum=sprintf('%.5d',num);
eta1=load([fdir 'Grd01_eta_' fnum]);
eta2=load([fdir 'Grd02_eta_' fnum]);
eta3=load([fdir 'Grd03_eta_' fnum]);
eta4=load([fdir 'Grd04_eta_' fnum]);
eta5=load([fdir 'Grd05_eta_' fnum]);

eta1(dep1<0)=NaN;
eta2(dep2<0)=NaN;
eta3(dep3<0)=NaN;
eta4(dep4<0)=NaN;
eta5(dep5<0)=NaN;

% subplot(length(nfile),1, num)

set(figure(1),'position',[50 50 620 500]) %figure default
pcolor(Xgrid1,Ygrid1,eta1),shading interp
colormap jet
hold on

pcolor(Xgrid2,Ygrid2,eta2),shading interp
pcolor(Xgrid3,Ygrid3,eta3),shading interp
pcolor(Xgrid4,Ygrid4,eta4),shading interp
pcolor(Xgrid5,Ygrid5,eta5),shading interp
% caxis([-0.1 0.1])
title([' Time = ', num2str(time(num+1)) ' hr '])

plot( [Xgrid2(1,1) Xgrid2(1,m2) Xgrid2(n2,m2) Xgrid2(n2,1) Xgrid2(1,1)], ...
    [Ygrid2(1,1) Ygrid2(1,m2) Ygrid2(n2,m2) Ygrid2(n2,1) Ygrid2(1,1)], ...
    'k--', 'LineWidth', 1.5 )

plot([Xgrid3(1,1) Xgrid3(1,m3) Xgrid3(n3,m3) Xgrid3(n3,1) Xgrid3(1,1)], ...
    [Ygrid3(1,1) Ygrid3(1,m3) Ygrid3(n3,m3) Ygrid3(n3,1) Ygrid3(1,1)], ...
    'k--', 'LineWidth', 1.5 )

plot([Xgrid4(1,1) Xgrid4(1,m4) Xgrid4(n4,m4) Xgrid4(n4,1) Xgrid4(1,1)], ...
    [Ygrid4(1,1) Ygrid4(1,m4) Ygrid4(n4,m4) Ygrid4(n4,1) Ygrid4(1,1)], ...
    'k--', 'LineWidth', 1.5 )

plot([Xgrid5(1,1) Xgrid5(1,m5) Xgrid5(n5,m5) Xgrid5(n5,1) Xgrid5(1,1)], ...
    [Ygrid5(1,1) Ygrid5(1,m5) Ygrid5(n5,m5) Ygrid5(n5,1) Ygrid5(1,1)], ...
    'k--', 'LineWidth', 1.5 )

caxis([-0.2 0.2])

% if num<=6
% caxis([-0.1 0.1])
% else
% caxis([-0.05 0.05])
% end

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')
cbar=colorbar;
set(get(cbar,'ylabel'),'String','\eta (m) ')
set(gcf,'Renderer','zbuffer');

axis equal; axis tight;

pname = ['surf_',num2str(num),'.png'];
set(gcf,'PaperPositionMode', 'auto');
print ('-dpng', pname);

end
%print -djpeg eta_inlet_shoal_irr.jpg