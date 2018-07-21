% This is to plot the profile of waveheights for test of berkhoffshoaling.

%------------------ dx=0.05;  dy=0.1
% section 1) x:=1m  	i=280
% section 2) x:=3m  	i=320
% section 3) x:=5m  	i=360
% section 4) x:=7m  	i=400
% section 5) x:=9m  	i=440
%------
% section 6) y:=-2m 	j=80
% section 7) y:= 0m 	j=100
% section 8) y:=+2m 	j=120
% center of the shoal : (x,y) = (260,100)
% ---------------------------------------

%------------------ dx=0.025;  dy=0.05
% section 1) x:=1m  	i=560
% section 2) x:=3m  	i=640
% section 3) x:=5m  	i=720
% section 4) x:=7m  	i=800
% section 5) x:=9m  	i=880
%------
% section 6) y:=-2m 	j=160
% section 7) y:= 0m 	j=200
% section 8) y:=+2m 	j=240
% center of the shoal : (x0,y0) = (520,200)
% ---------------------------------------

clear all;  close all;

%-----------------------
mx=1200;  ny=400;
dx=0.025;  dy=0.05;
x0=520;     y0=200;    % center of the shoal
%-----------------------

x(1:mx)=0.0;
y(1:ny)=0.0;

xshift=0;
yshift=2;

for i=1:mx
    x(i)=(i-1-x0)*dx;
end
for j=1:ny
    y(j)=(j-1-y0)*dy;
end

height=load('height_full.dat');      %computated wave height 
%load depth.dat;   %depth file
%load tgzeta7.dat; %a surface elevation file
wh = height;

load section1.dat;    %experimental data with amplitude in mm
load section2.dat;
load section3.dat;
load section4.dat;
load section5.dat;
load section678.dat; 

%set(figure(1),'position',[50 50 600 600])
figure(1)
% section 1----
  subplot(4,2,1);
  plot(section1(:,1),section1(:,2)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(y,wh(560,:)*1000, 'LineWidth', 1.2);
%   plot(y,wh(560+xshift,:)*1000);
  %plot(y,wh(540+xshift,:)*1000);
  axis([-5,5,10,80]);
  xlabel('y (m)');
  ylabel('H(mm)')

% section 2----
  subplot(4,2,2);
  plot(section2(:,1),section2(:,2)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(y,wh(640,:)*1000, 'LineWidth', 1.2);
  axis([-5,5,10,80]);
  xlabel('y (m)');
  ylabel('H(mm)')

% section 3----
  subplot(4,2,3);
  plot(section3(:,1),section3(:,2)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(y,wh(720,:)*1000, 'LineWidth', 1.2);
  axis([-5,5,10,120]);
  xlabel('y (m)');
  ylabel('H(mm)')

% section 4----
  subplot(4,2,4);
  plot(section4(:,1),section4(:,2)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(y,wh(800,:)*1000, 'LineWidth', 1.2);
  axis([-5,5,-5,120]);
  xlabel('y (m)');
  ylabel('H(mm)')

% section 5----
  subplot(4,2,5);
  plot(section5(:,1),section5(:,2)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(y,wh(880,:)*1000, 'LineWidth', 1.2);
  axis([-5,5,5,100]);
  xlabel('y (m)');
  ylabel('H(mm)')

% section 6----
  subplot(4,2,6);
  plot(-section678(:,1),section678(:,2)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(x,wh(:,160+yshift)*1000, 'LineWidth', 1.2);
  axis([0,11,5,120]);
  xlabel('x (m)');
  ylabel('H(mm)')

% section 7----
  subplot(4,2,7);
  plot(-section678(:,1),section678(:,3)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(x,wh(:,200+yshift)*1000, 'LineWidth', 1.2);
  axis([0,11,5,120]);
  xlabel('x (m)');
  ylabel('H(mm)')


% section 8----
  subplot(4,2,8);
  plot(-section678(:,1),section678(:,4)*2,'ko','MarkerSize',3.5, 'LineWidth', 1.0);
  hold;
  plot(x,wh(:,240+yshift)*1000, 'LineWidth', 1.2);
  axis([0,11,5,120]);
  xlabel('x (m)');
  ylabel('H(mm)')

  outfig=['.\result.eps'];
  set(gcf, 'PaperPositionMode', 'auto');
  print('-depsc2',outfig)

% The timegage 7 and topograph
%  figure(2)
%  pcolor(tgzeta7)
%  shading interp;
%  hold;
%  contour(-depth/5,10)
%
