clear all; clc;

%==============
dx1=0.2;     dy1=0.2;
%==============
Nghost=4;

fdir='/Users/fengyanshi15/tmp2/';
fdir1='../';

% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

grid=load(['../04_serial_VBtest2/subgrid_info.txt']);
%dep=load([fdir 'Grd01_dep.out']);

n1=126;
m1=165;
%[n1, m1]=size(dep);
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

%-------------------------------------------------------------

%time=load([fdir 'Grd01_track.txt']);
%nstart=input('nstart=');
%nend=input('nend=');

nstart=2;
nend=150;

icount=0;
for num=nstart:1:nend
icount=icount+1;

fname=sprintf('%.5d',num);
eta1=load([fdir 'eta_' fname]);

time(num)=num*0.1;

clf
pcolor(Xgrid1,Ygrid1,eta1),shading interp
hold on

xlabel('x (m) ', 'fontsize', 12)
ylabel('y (m) ', 'fontsize', 12)

colormap jet;
colorbar
caxis([-0.02,0.02])
axis([10 27.5 0 25])


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
