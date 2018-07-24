clear all; clc;

%==============
dx1=0.2;     dy1=0.2;
%==============
Nghost=4;

fdir='/Users/fengyanshi15/tmp3/';
fdir1='../';

% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

grid=load(['../04_serial_VBtest2/subgrid_info.txt']);
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
x2 = x2 + (trackx2(1)-1)*dx1 - (Nghost*dx2);
y2 = y2 + (tracky2(1)-1)*dy1 - (Nghost*dy2);
[Xgrid2, Ygrid2]=meshgrid(x2,y2);
%-------------------------------------------------------------

m3=grid(2,1);
n3=grid(2,2);
isk3=grid(2,3);
trackx3=grid(2,4);
tracky3=grid(2,5);

dx3=dx2/isk3;
dy3=dy2/isk3;

x3 = 0 : dx3 : (m3-1)*dx3;
y3 = 0 : dy3 : (n3-1)*dy3;

x3 = x3 + min( x2 ) + (trackx3-1)*dx2 - (Nghost*dx3);
y3 = y3 + min( y2 ) + (tracky3-1)*dy2 - (Nghost*dy3);
[Xgrid3, Ygrid3]=meshgrid(x3,y3);
%-------------------------------------------------------------

time=load([fdir 'Grd01_track.txt']);
%nstart=input('nstart=');
%nend=input('nend=');

nstart=2;
nend=150;

icount=0;
for num=nstart:1:nend
icount=icount+1;

fname=sprintf('%.5d',num);
eta1=load([fdir 'Grd01_eta_' fname]);
eta2=load([fdir 'Grd02_eta_' fname]);
eta3=load([fdir 'Grd03_eta_' fname]);

clf
pcolor(Xgrid1,Ygrid1,eta1),shading interp
hold on

xlabel('x (m) ', 'fontsize', 12)
ylabel('y (m) ', 'fontsize', 12)

pcolor(Xgrid2,Ygrid2,eta2),shading interp
pcolor(Xgrid3,Ygrid3,eta3),shading interp
%caxis([-0.1 0.1])
colormap jet;
colorbar
caxis([-0.02,0.02])
axis([10 27.5 0 25])

plot( [Xgrid2(1,1) Xgrid2(1,m2) Xgrid2(n2,m2) Xgrid2(n2,1) Xgrid2(1,1)], ...
    [Ygrid2(1,1) Ygrid2(1,m2) Ygrid2(n2,m2) Ygrid2(n2,1) Ygrid2(1,1)], ...
    'k--', 'LineWidth', 1.7 )

plot([Xgrid3(1,1) Xgrid3(1,m3) Xgrid3(n3,m3) Xgrid3(n3,1) Xgrid3(1,1)], ...
    [Ygrid3(1,1) Ygrid3(1,m3) Ygrid3(n3,m3) Ygrid3(n3,1) Ygrid3(1,1)], ...
    'k--', 'LineWidth', 1.3)

M(:,icount)=getframe(gcf);

title (['Time = ',num2str(time(num), '%5.2f'), ' (s)'],'fontsize',11);
box on;

pause(0.1)

    currframe=getframe(gcf);
    writeVideo(vidObj,currframe);  % Get each recorded frame and write it to filename defined above


pname = ['surf_',num2str(fname),'.png'];
%print ('-dpng', pname);

% pname = ['surf_',num2str(fname),'.bmp'];
% print ('-dbmp','-r300', pname);

end
close(vidObj)
