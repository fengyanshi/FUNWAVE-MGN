clear all; clc;

%==============
dx1=10;     dy1=10;
%==============
Nghost=4;

fdir='/Users/fengyanshi15/tmp1/';
fdir1='../';

grid=load(['subgrid_info.txt']);
dep=load([fdir 'Grd01_dep.out']);

[n1, m1]=size(dep);
x1 = 0 : dx1 : (m1-1)*dx1;
y1 = 0 : dy1 : (n1-1)*dy1;
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

% x2 = x1( trackx2(1) ) + x2;
% y2 = y1( tracky2(1) ) + y2;
x2 = x2 + (trackx2(1)-1)*dx1 - (3*Nghost*dx2);
y2 = y2 + (tracky2(1)-1)*dy1 - (3*Nghost*dy2);
[Xgrid2, Ygrid2]=meshgrid(x2,y2);
%-------------------------------------------------------------

time=load([fdir 'Grd01_track.txt']);
nstart=input('nstart=');
nend=input('nend=');

icount=0;
for num=nstart:1:nend
icount=icount+1;

fname=sprintf('%.5d',num);
eta1=load([fdir 'Grd01_eta_' fname]);
eta2=load([fdir 'Grd02_eta_' fname]);

clf
pcolor(Xgrid1,Ygrid1,eta1),shading interp

hold on

%pcolor(Xgrid2,Ygrid2,eta2),shading interp
%caxis([-0.1 0.1])
colormap jet;
colorbar

plot( [Xgrid2(1,1) Xgrid2(1,m2) Xgrid2(n2,m2) Xgrid2(n2,1) Xgrid2(1,1)], ...
    [Ygrid2(1,1) Ygrid2(1,m2) Ygrid2(n2,m2) Ygrid2(n2,1) Ygrid2(1,1)], ...
    'k--', 'LineWidth', 1.7 )
M(:,icount)=getframe(gcf);

title (['Time = ',num2str(time(num), '%5.2f'), ' (s)'],'fontsize',11);
box on;

pause(0.5)

pname = ['surf_',num2str(fname),'.png'];
%print ('-dpng', pname);

% pname = ['surf_',num2str(fname),'.bmp'];
% print ('-dbmp','-r300', pname);

end