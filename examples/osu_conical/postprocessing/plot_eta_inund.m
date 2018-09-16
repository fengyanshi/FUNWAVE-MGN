%clear all; clc; close all; fclose all; 

fdir = '/Users/fengyanshi15/tmp1/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depthname=[fdir 'Grd01_dep.out'];
depthname2=[fdir 'Grd02_dep.out'];
depthname3=[fdir 'Grd03_dep.out'];
grid=load(['subgrid_info.txt']);

mindep=0.005;
mx=125;        ny=66;
dx=0.4;          dy=0.4;
Nghost=4;
x1=0;         x2=(mx-1)*dx;
y1=0;         y2=(ny-1)*dy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x1 : dx : x2;
y = y1 : dy : y2;
x = x - 5;
y = y - 13;
[xgrid,ygrid]=meshgrid(x,y);  clear x1 y1 x2 y2;
%%%%%%%%%

m2=grid(1,1);
n2=grid(1,2);
isk2=grid(1,3);
trackx2=grid(1,4);
tracky2=grid(1,5);

dx2=dx/isk2;
dy2=dy/isk2;

%-------------------------------------------------------
x2 = 0 : dx2 : (m2-1)*dx2;
y2 = 0 : dy2 : (n2-1)*dy2;
x2 = x2 + (trackx2-Nghost-1)*dx + Nghost*dx2;
y2 = y2 + (tracky2-Nghost-1)*dy + Nghost*dy2;
%-------------------------------------------------------

x2 = x2 - 5;
y2 = y2 - 13;
[xgrid2, ygrid2]=meshgrid(x2,y2);
%%%%%%%%%

m3=grid(2,1);
n3=grid(2,2);
isk3=grid(2,3);
trackx3=grid(2,4);
tracky3=grid(2,5);

dx3=dx2/isk3;
dy3=dy2/isk3;

%-------------------------------------------------------
x3 = 0 : dx3 : (m3-1)*dx3;
y3 = 0 : dy3 : (n3-1)*dy3;
x3 = x3 + min( x2 ) + (trackx3-Nghost-1)*dx2 + Nghost*dx3;
y3 = y3 + min( y2 ) + (tracky3-Nghost-1)*dy2 + Nghost*dy3;
%-------------------------------------------------------

[xgrid3, ygrid3]=meshgrid(x3,y3);
%%%%%%%%%

time=load([fdir 'Grd01_track.txt']);
tm=time(:,1);  clear time

dep=load( depthname );
dep2=load( depthname2 );
dep3=load( depthname3 );

num_file=input('file number');
%for kk=1:1:length(tm);
for kk=num_file:num_file;
    fnum=kk;
    
    fname = ['Grd01_eta_',num2str(fnum,'%05d')];
    eta=load( [ fdir, fname ] );
    
    ind=find( eta + dep < mindep );
    eta(ind) = NaN;
    %%%%%
    
    fname = ['Grd02_eta_',num2str(fnum,'%05d')];
    eta2=load( [ fdir, fname ] );
    
    ind=find( eta2 + dep2 < mindep );
    eta2(ind) = NaN;    
    %%%%%
    
    fname = ['Grd03_eta_',num2str(fnum,'%05d')];
    eta3=load( [ fdir, fname ] );
    
    ind=find( eta3 + dep3 < mindep );
    eta3(ind) = NaN;        
    %%%%%
    
    %figure(kk)
    clf

     surface(xgrid,ygrid,-dep,...
         'FaceColor',[0.68 0.46 0.0],'EdgeColor','none','CDataMapping','direct',...
         'DiffuseStrength',0.5)

     hold on        

     surface(xgrid2,ygrid2,-dep2,...
         'FaceColor',[0.68 0.46 0.0],'EdgeColor','none','CDataMapping','direct',...
         'DiffuseStrength',0.5)

     surface(xgrid3,ygrid3,-dep3,...
         'FaceColor',[0.68 0.46 0.0],'EdgeColor','none','CDataMapping','direct',...
         'DiffuseStrength',0.5)

     surface(xgrid,ygrid,eta,...
         'FaceColor',[0.2 0.5 1],'EdgeColor','none','CDataMapping','direct',...
         'DiffuseStrength',0.5)
     lightangle(-30,30)

     surface(xgrid2,ygrid2,eta2,...
         'FaceColor',[0.2 0.5 0.8],'EdgeColor','none','CDataMapping','direct',...
         'DiffuseStrength',0.5)
     lightangle(-30,30)
    
     surface(xgrid3,ygrid3,eta3,...
         'FaceColor',[0.2 0.5 0.6],'EdgeColor','none','CDataMapping','direct',...
         'DiffuseStrength',0.5)
     lightangle(-30,30)
    
     plot( [ xgrid2(1,1) xgrid2(1,m2) xgrid2(n2,m2) xgrid2(n2,1) xgrid2(1,1) ], ...
         [ ygrid2(1,1) ygrid2(1,m2) ygrid2(n2,m2) ygrid2(n2,1) ygrid2(1,1)], ...
         'w--', 'LineWidth', 1.5 )    

     plot( [ xgrid3(1,1) xgrid3(1,m3) xgrid3(n3,m3) xgrid3(n3,1) xgrid3(1,1) ], ...
         [ ygrid3(1,1) ygrid3(1,m3) ygrid3(n3,m3) ygrid3(n3,1) ygrid3(1,1)], ...
         'w--', 'LineWidth', 1.5 )
     
     %time = tm( kk );
     %title (['time = ',num2str( time, '%5.2f' ), ' (s)'],'fontsize',11);
    
%     view(-3,34)    
     view(30,58)    
    
    xlim([min(xgrid(:)) max(xgrid(:))]);
    ylim([min(ygrid(:)) max(ygrid(:))]);
    zlim([-1.0 0.5])
        
    xlabel('x (m)','fontsize',12)
    ylabel('y (m)','fontsize',12)
    zlabel('z (m)','fontsize',12)
    
    set(gca, 'fontsize', 10);
    
    pname = ['surfeta_', num2str( kk, '%5.5i' ) ,'.png'];
    %print ('-dpng', pname);
    %hold off; close all; fclose all; 
end