% $Header: /u/gcmpack/MITgcm/verification/advect_xz/input/gendata.m,v 2.2 2011/12/05 22:03:52 jmc Exp $
% $Name: checkpoint65r $
clear all;
addpath(genpath('/glade/u/home/liuchihl/TOOLS'));

% This is a matlab script that generates the input data
prec='real*8';
ieee='b';
w2file=1; %- do write to file (1=y/0=no)

%- Dimensions of grid
%nx=1377;
nx=1200;
ny=1;
nz=100;

%- Horizontal & vertical resolution (m):

%dx = [340:-2:20,20*ones(1,1075),20:2:300]; %dy=dx; 
dx = [300:-2:20,20*ones(1,918),20:2:300];
dy=20;
dz=[ones(1,30)*3,ones(1,50)*3,4:23];

fid=fopen('delX','w','b');fwrite(fid,dx,'real*8');fclose(fid);
fid=fopen('delZ','w','b');fwrite(fid,dz,'real*8');fclose(fid);

% full size of the domain:
Lx=sum(dx); % x domain
H=sum(dz) ; % vertical depth

%- grid point coordinate :
% f indicates faces; c indicates cell center
xf(1)=0;
for i=1:nx
xf(i+1)=xf(i)+dx(i);
end
xc=(xf(1:nx)+xf(2:nx+1))/2;

yf(1)=0;
for i=1:ny
yf(i+1)=yf(i)+dy;
end
yc=(yf(1:ny)+yf(2:ny+1))/2;

zf(1)=0;
for i=1:nz
    zf(i+1)=zf(i)+dz(i);
end
zc=(zf(1:nz)+zf(2:nz+1))/2;



%- bathymetry :

% We use a simple Gaussian Bump as an example

dh = 0.8235*H;              % height of bump
 L = Lx/25;                  % scale of the bump
x0 = mean(xc); 
zb = -H+dh*exp( -0.9*( ((xc-x0).^2)/(L^2) ) );
% close all;figure; plot(xc,zb)
zb3D = zb*ones(1,ny);       %3D uniform bathymetry

% save bathymetry into binary file
wf=' write file: ';
if w2file;
 bname='bathy_slope.bin';
 fid=fopen(bname,'w',ieee); fwrite(fid,zb3D,prec); fclose(fid); 
 fprintf([wf,bname,'\n']);
end

save topo_sill zb zc xc yc

%- partial cell filling (hFac):
% hFac=ones(nx,1)*zf(1,[2:nz+1]);
% hFac=max(hFac,zb*ones(1,nz));           
% hFac=ones(nx,1)*zf(1,[1:nz])-hFac;
% hFac=hFac/dz; hFac=max(0,hFac);
% 
% rhFc=reshape(hFac,[nx*nz 1]); 
% [IK]=find(rhFc > 0); rhFc(IK)=1./rhFc(IK);
% rhFc=reshape(rhFc,[nx nz]);
% rhFc(rhFc>1)=1;
% for i=1:ny
% PC(:,i,:) = rhFc;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- velocity:

%%% Input velocity at the boundary
ua = 0.5*(tanh((-zc+100)/120)+1);
dt=1/60;        % 1 min in hour
t=dt:dt:72;     % give a 3 day tidal current
% omega = 1.41*1e-4; %M2 tide frequency
T = 12.41667;   % tidal period in hours
U0 = .12;       % tidal amplitude   

bt_timedelay = Lx/sqrt(9.81*-mean(zb3D(:,1))); % barotropic wave delay
    for j=1:ny
        for tt=1:length(t)
% inflow velocity: ua+tide (starting from the low tide)
Uvel_in(1,j,:,tt) = ua+U0*sin(2*pi/T*(t(tt)-t(1))-pi/2 );   
Uvel_out(1,j,:,tt) = ua+U0*sin(2*pi/T*(t(tt)-t(1)-bt_timedelay/3600)-pi/2);   % current velocity + tidal velocity
        end
    end
% spin-up function f(t)
stime = 1381;                      % spinup time 1 day
f = ones(1,ny,nz,length(t));
f(1,ny,:,1:stime) = ones(nz,1)*linspace(0,1,stime); % t(313) is at the half tidal cycle
Uvel_in = Uvel_in.*f; Uvel_in(1,:,:,end) = 0;
Uvel_out = Uvel_out.*f; Uvel_out(1,:,:,end) = 0;
figure;
til = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(squeeze(Uvel_in(1,1,:,stime)),-zc,'-','linewidth',2);
hold on;plot(squeeze(Uvel_out(1,1,:,stime)),-zc,'-','linewidth',2);
xlabel('U(z) m/s'); ylabel('Depth (m)');

nexttile

plot(t,squeeze(Uvel_in(1,1,1,:)),'-','linewidth',2);
hold on; plot(t,squeeze(Uvel_out(1,1,1,:)),'-','linewidth',2)
xlabel('t (hr)'); ylabel('U_{obcs} (m/s)');
% hold on; plot(squeeze(Uvel(1,:,[45])),zc);

%save ua to binary (Uvel is the boundary velocity)
uname = 'Uvel_in.bin';
fid=fopen(uname,'w',ieee); fwrite(fid,Uvel_in,prec); fclose(fid);
fprintf([wf,uname,'\n']);
uname = 'Uvel_out.bin';
fid=fopen(uname,'w',ieee); fwrite(fid,Uvel_out,prec); fclose(fid);
fprintf([wf,uname,'\n']);

%%% initial velocity 
[~,crest] = max(zb3D(:,1));      % the position of the crest

Uini2D = ones(nx,1)*squeeze(Uvel_in(1,1,:,1))';                   % x-z velocity field without adding tide

for k=1:nz
    for j=1:ny
        for i=1:nx
       Uini(i,j,k,1) = Uini2D(i,k);       %3D initial velocity
        end
    end
end
% initial random noise in the y direction
% Uini = Uini+(rand(nx,ny,nz)-.5)*.01;
% Uini = Uini; 
% 
% Vini = (rand(nx,ny,nz)-.5)*.01;		

Uininame = 'Uini.bin';
fid=fopen(Uininame,'w',ieee); fwrite(fid,Uini,prec); fclose(fid);
fprintf([wf,Uininame,'\n']);
% Vininame = 'Vini.bin';
% fid=fopen(Vininame,'w',ieee); fwrite(fid,Vini,prec); fclose(fid);
% fprintf([wf,Vininame,'\n']);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here we set the temperature, salinity and other tracers below (if any)

%- stratification
g=9.81;
talpha=2.0e-4;
N2 = 1.44e-4;
Tz = N2/(g*talpha);
%- temperature
% temperature at the boundary
% Tref = Tz*(-zc)-mean(Tz*(-zc));
% Tref = Tz*(-zc)+28;
Tref = linspace(28,10,nz);
% xlim([5,30])
% hold on; plot(Tref,-zc);

% initial temperature in the 2D (x-z) domain
tref = ones(nx,1)*Tref;
tref(tref>max(Tref))=0;

%%%- salinity
% salinity at the boundary
Sref = 35*ones(1,nz);

% initial salinity in the 2D (x-z) domain
sref = ones(nx,1)*Sref;
sref = sref; sref(sref>max(Sref))=0;

% expand the T/S initial condition to a 3D matrix
clear Tini Sini
for k=1:nz
    for j=1:ny
        for i=1:nx
       Tini(i,j,k,1) = tref(i,k);       %3D initial T
       Sini(i,j,k,1) = sref(i,k);       %3D initial S
        end
    end
end

%- make boundary T,S as a 3D matrix (add time to fit tidal effect)
clear TT SS
for j=1:ny
    for tt=1:length(t)
        TT(1,j,:,tt) = Tref;
        SS(1,j,:,tt) = Sref;
    end
end
% save initial T,S as binary file
if w2file,
 tname='Tini_G.bin';sname='Sini_G.bin';
 fid=fopen(tname,'w',ieee); fwrite(fid,Tini,prec); fclose(fid); 
 fid=fopen(sname,'w',ieee); fwrite(fid,Sini,prec); fclose(fid); 
 fprintf([wf,tname,'\n']);fprintf([wf,sname,'\n']);
end
% save the T/S value for open boundary condition 
tOBname='OB_T.bin';%boundary temperature
fid=fopen(tOBname,'w',ieee); fwrite(fid,TT,prec); fclose(fid); 
fprintf([wf,tOBname,'\n']);
sOBname='OB_S.bin';%boundary salinity
fid=fopen(sOBname,'w',ieee); fwrite(fid,SS,prec); fclose(fid); 
fprintf([wf,sOBname,'\n']);


%%
fid = fopen('Uvel_in.bin')
Uvel = fread(fid, 'real*8', 'ieee-be');
fclose(fid);
plot(Uvel)
%% - make plots to check:
close all;
% figure('position',[100 100 800 600]);clf;
% pcolor(xf(1:nx)*1e-3,zc,squeeze(Uini(:,1,:))'); shading flat;colorbar;
% hold on ; plot(xc*1e-3,zb,'r-'); hold off ; grid ; ylim([-H 0]);
% title('U velocity [m/s]');

figure('position',[100 100 800 600]);clf;

til = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
 pcolor(xc*1e-3,-zc,tref'); shading flat;colorbar;caxis([min(Tref) max(Tref)]);
 hold on ; plot(xc*1e-3,zb,'r-'); hold off ; grid; ylim([-H 0]);
 title('Initial Theta');
nexttile
pcolor(xf(1:nx)*1e-3,-zc,squeeze(Uini(:,1,:))'); shading flat;colorbar;
 hold on ; plot(xc*1e-3,zb,'r-'); hold off ; grid; ylim([-H 0]);
title('U velocity [m/s]');

return
