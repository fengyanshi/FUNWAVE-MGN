clear all
fdir='/Users/fengyanshi15/tmp1/';
Nghost=4;
totalres=48;
timelag=3600;
dep=load('../external_files/depth_30min.txt');

grid=load(['subgrid_info.txt']);

[n1,m1]=size(dep);
dx1=0.5;
dy1=0.5;
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

x2 = x2 + (trackx2-1)*dx1 - (Nghost*dx2);
y2 = y2 + (tracky2-1)*dy1 - (Nghost*dy2);
[Xgrid2, Ygrid2]=meshgrid(x2,y2);
%-------------------------------------------------------------

% nfile=[2 9];
% hr={'1' '8'};

% wid=5;
% len=7;
% set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
% clf

time=load([fdir 'Grd01_track.txt']);
time=time+timelag;
time=time/3600;

for num=1:totalres;
    
fnum=sprintf('%.5d',num);
eta1=load([fdir 'Grd01_eta_' fnum]);
eta2=load([fdir 'Grd02_eta_' fnum]);

eta1(dep<0)=NaN;

% subplot(length(nfile),1, num)

set(figure(1),'position',[50 50 620 500]) %figure default
pcolor(Xgrid1,Ygrid1,eta1),shading interp
colormap jet
hold on

pcolor(Xgrid2,Ygrid2,eta2),shading interp
% caxis([-0.1 0.1])
%title([' Time = ', num2str(time(num)) ' hr '])

caxis([-0.5 0.5])

plot( [Xgrid2(1,1) Xgrid2(1,m2) Xgrid2(n2,m2) Xgrid2(n2,1) Xgrid2(1,1)], ...
    [Ygrid2(1,1) Ygrid2(1,m2) Ygrid2(n2,m2) Ygrid2(n2,1) Ygrid2(1,1)], ...
    'k--', 'LineWidth', 1.5 )

if num<=7
caxis([-0.1 0.1])
else
caxis([-0.03 0.03])
end

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