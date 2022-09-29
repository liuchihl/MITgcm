clear all;load uvt2.mat;
addpath(genpath('/home1/chihlun/MITgcm_c65r/mhlib'));
load('/home1/chihlun/Data_base/200mData/topo.mat');
nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.15,122.418,nx);
y_new = linspace(24.29,24.9,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';

close all;
dim=size(U);
U_z=nan(dim(1),dim(2),dim(3),dim(4));
% x_range = 790:970;
x_range = 700:970;

for i=[1:length(itv)]
for j=x_range %X-axis
ut=squeeze(U(j,:,:,i));
U_z(j,:,:,i)=mmderiv(-zc,ut')';
end
end

%- plot shear, V, W altogether 
%  U_z(U_z==0)=nan;T(T==0)=nan;U(U==0)=nan;W(W==0)=nan;V(V==0)=nan;
%  set(gca,'nextplot','replacechildren');
 temp = VideoWriter('Av.mp4');
 temp.FrameRate = 5;
 temp.Quality = 100;
 open(temp);

% itv = ITS;
% for i=1:length(itv)
% for i=[1:382]
% for i=[1:length(itv)]
    for i=1:length(itv)
figure;
%%%%%%%%%%%% figure size
% set(gcf,'units','centimeters','paperunits','centimeters')
% set(gcf,'PaperType','A4');
% % pp=[0.63 0.9 27.5 28];
% pp=[0.63 0.9 31 28];
% ps=[0 0 pp(3)/1.1 pp(4)/1.1];
% set(gcf,'paperposition',pp)
% set(gcf,'position',ps)
figset(2,1);
%%%%%%%%%%%%

%- T
ax1 = axes('position',[0.1 0.68 .8 .3]);
[c,hh] = contourf(xc(x_range)*1e-3,zc,squeeze(KLviscAr(x_range,1,:,i))',[1e-5:0.001:.9,1:0.1:3.5]);
set(hh,'edgecolor','none');
caxis([1e-4 5e-1]);
% colormap(ax1,(flipud(cbrewer('seq','YlGnBu',100))));
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[15:.8:28],'color','w');
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[15:.3:28],'color',[40 40 40]/225);
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[22 22],'color','y');
set(hh,'linewidth',2);
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[17 17],'color','g');
set(hh,'linewidth',2);

z1 = colorbar('location','eastoutside');
set(z1,'position',[0.92 0.76 .02 0.22]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
% set(gca,'xticklabel',[]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
% colormap(ax1,(cbrewer('div','RdYlBu',100)));
cmocean('curl')
% colormap(ax1,'default');
% colormap(ax1,flipud(brewermap([],'Spectral')))
% load MODIS_colorbar.mat;MODIS = colormap;clear colormap;
% colormap(ax1,(MODIS));
% xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);

% set(gca,'xtick',[36.5:2:48.48]);
% set(gca,'xticklabel',{'0','2','4','6','8','10','12'});
set(gca,'fontsize',15);
tit = title(sprintf('Time=%.1f minutes',2*i));
set(tit,'fontsize',10)
text(41,-220,'T (^oC)','fontsize',15,'fontname','times');

%- U
ax2 = axes('position',[0.1 0.522 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((U(x_range,1,:,i)))',[-1:0.01:2]);
caxis([-.5 1.5]);
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((U((x_range),1,:,i)))',[-1:0.2:2],'k');
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((U((x_range),1,:,i)))',[0 0],'r');
set(hh,'linewidth',2);
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
set(gca,'xticklabel',[]);
myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
set(gca,'tickdir','out')
colormap(ax2,flipud(cbrewer('div','RdGy',100)));
z2 = colorbar('location','eastoutside');
set(z2,'position',[0.92 0.522 .02 0.22]);
hold on;fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'fontsize',15)
text(41,-220,'U (ms^-^1)','fontsize',15,'fontname','times');
%- W
ax3 = axes('position',[0.1 0.284 .8 .22]);
 [c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((W(x_range,1,:,i)))',[-1:0.005:1]);
 caxis([-.5 .5]);
 set(hh,'edgecolor','none');
 z3 = colorbar('location','eastoutside');
 set(z3,'position',[0.92 0.284 .02 0.22]);
set(gca,'tickdir','out')
colormap(ax3,flipud(cbrewer('div','RdGy',100)));
% hold on;
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((W((x_range),1,:,i)))',[-1:0.5:2],'k');
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
set(gca,'xticklabel',[]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
hold on; plot(xc(815)*1e-3,0,'or','markersize',15);
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'fontsize',15)
text(41,-220,'W (ms^-^1)','fontsize',15,'fontname','times');

%- shear square
ax4 = axes('position',[0.1 0.046 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
    [0:0.01:1,1:3:40]);
set(hh,'linewidth',0.1)
caxis([0 1]);%1
set(hh,'edgecolor','none');
set(gca,'xtick',[0:5:60]);
set(gca,'tickdir','out')
% cmocean('delta');
cmocean('-deep');
% colormap(ax4,'default');
% load MODIS_colorbar.mat;MODIS = colormap;clear colormap;
% colormap(ax4,(MODIS));
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
    [-0.2:0.1:0.2],'color','k');
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);

% set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
% set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
set(gca,'xticklabel',[{0} {0.5} {1.0} {1.5} {2.0} {2.5} {3.0} {3.5} {4.0} {4.5} {5.0} ]);
myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
z4 = colorbar('location','eastoutside');
set(z4,'position',[0.92 0.046 .02 0.22]); 
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
ylabel('Depth (m)','fontsize',15,'fontname','times')
xxx = xlabel('X (km)','fontsize',15,'fontname','times');
% get(xxx);
set(xxx,'position',get(xxx,'position')+[2.8 22 1.0]);
% set(xxx,'position',[0.0 .5 0.5]);

set(gca,'fontsize',15);
set(gcf,'color','w');
text(41,-220,'S^2 (10^-^3s^-^2)','fontsize',15,'fontname','times');

ax5 = axes('position',[0.66 0.33 .226 .164]);
U_interest = squeeze(U(815,1,1,:));
U_lowpass = lpass(U_interest,2/60,2,2);
% figure;
% plot((itv-itv(1))/3600,U_lowpass);
% hold on;plot((itv-itv(1))/3600,U_interest)
plot((itv-itv(1))/3600,U_lowpass,'k-','linewidth',2);
hold on; plot((itv(i)-itv(1))/3600,U_lowpass(i),'ro','markersize',15);
axis([0 25 0.5 2]);
ylabel('U(m/s)');xlabel('Time(hour)');grid on;
set(gca,'fontsize',10);
%------------------ time lag is too large for the boundary U
t=linspace(0,24,length(itv));
Period = 12.41667; %tidal period in hours
ua=(-0.15*atan((-zc-200)/100)+0.5)*.8;% Along stream velocity
for tt=1:length(t)
Uvel(1,:,tt) = ua+.5*ua.*sin(2*pi/Period*t(tt));
end
plot(t,squeeze(Uvel(1,1,:)),'linewidth',2,'color','k');
hold on; plot(t(i),squeeze(Uvel(1,1,i)),'ro');
set(gcf,'color','w');
xlabel('Time (hour)');ylabel('U (m/s)'); set(gca,'fontsize',15);grid on;
%---------------------
drawnow 
 frame = getframe(gcf);
 if size(frame.cdata,1)~=926||size(frame.cdata,2)~=1065
%      figset(2,1);
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 31 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp)
set(gcf,'position',ps)
 frame = getframe(gcf);

 end
 
 writeVideo(temp,frame);
 close;
end
close(temp);



%%

clear all;load uvt2.mat;
addpath(genpath('/home1/chihlun/MITgcm_c65r/mhlib'));
load('/home1/chihlun/Data_base/200mData/topo.mat');
nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.15,122.418,nx);
y_new = linspace(24.29,24.9,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';

close all;
dim=size(U);
U_z=nan(dim(1),dim(2),dim(3),dim(4));
% x_range = 790:970;
x_range = 700:970;

for i=[1:length(itv)]
for j=x_range %X-axis
ut=squeeze(U(j,:,:,i));
U_z(j,:,:,i)=mmderiv(-zc,ut')';
end
end

%- plot shear, V, W altogether 
%  U_z(U_z==0)=nan;T(T==0)=nan;U(U==0)=nan;W(W==0)=nan;V(V==0)=nan;
%  set(gca,'nextplot','replacechildren');
 temp = VideoWriter('TUW_hd_time.avi');
 temp.FrameRate = 2;
 temp.Quality = 100;
 open(temp);

% itv = ITS;
% for i=1:length(itv)
% for i=[1:382]
% for i=[1:length(itv)]
    for i=1:length(itv)
figure;
%%%%%%%%%%%% figure size
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
% pp=[0.63 0.9 27.5 28];
pp=[0.63 0.9 31 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp)
set(gcf,'position',ps)
%%%%%%%%%%%%

%- T
ax1 = axes('position',[0.1 0.76 .8 .22]);
% [c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[14:.1:28.5]);%caxi
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[14:.1:18.2,18.3:0.005:18.7,18.8:0.1:28.5]);%caxi
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[17:.2:22]);%caxi
% colormap(ax1,(flipud(cbrewer('seq','YlGnBu',100))));
% caxis([16 22]);
% set(gca,'color','k');
% hold on;
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[22:.4:28],'color',[40 40 40]/225 );
% hold on;
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[14:.4:16.9],'color',[40 40 40]/225);
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[15:.3:28],'color',[40 40 40]/225);
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[22 22],'color','y');
set(hh,'linewidth',2);
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[17 17],'color','g');
set(hh,'linewidth',2);
caxis([15 28]);

z1 = colorbar('location','eastoutside');
set(z1,'position',[0.92 0.76 .02 0.22]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
set(gca,'xticklabel',[]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
colormap(ax1,flipud(cbrewer('div','RdYlBu',100)));
% colormap(ax1,'default');
% colormap(ax1,flipud(brewermap([],'Spectral')))
% load MODIS_colorbar.mat;MODIS = colormap;clear colormap;
% colormap(ax1,(MODIS));
% xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);

% set(gca,'xtick',[36.5:2:48.48]);
% set(gca,'xticklabel',{'0','2','4','6','8','10','12'});
set(gca,'fontsize',15);
tit = title(sprintf('Time=%.1f minutes',2*i));
set(tit,'fontsize',10)
text(41,-220,'T (^oC)','fontsize',15,'fontname','times');

%- U
ax2 = axes('position',[0.1 0.522 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((U(x_range,1,:,i)))',[-1:0.01:2]);
caxis([-.5 1.5]);
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((U((x_range),1,:,i)))',[-1:0.2:2],'k');
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((U((x_range),1,:,i)))',[0 0],'r');
set(hh,'linewidth',2);
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
set(gca,'xticklabel',[]);
myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
set(gca,'tickdir','out')
colormap(ax2,flipud(cbrewer('div','RdGy',100)));
z2 = colorbar('location','eastoutside');
set(z2,'position',[0.92 0.522 .02 0.22]);
hold on;fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'fontsize',15)
text(41,-220,'U (ms^-^1)','fontsize',15,'fontname','times');
%- W
ax3 = axes('position',[0.1 0.284 .8 .22]);
 [c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((W(x_range,1,:,i)))',[-1:0.005:1]);
 caxis([-.5 .5]);
 set(hh,'edgecolor','none');
 z3 = colorbar('location','eastoutside');
 set(z3,'position',[0.92 0.284 .02 0.22]);
set(gca,'tickdir','out')
colormap(ax3,flipud(cbrewer('div','RdGy',100)));
% hold on;
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((W((x_range),1,:,i)))',[-1:0.5:2],'k');
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
% set(gca,'xticklabel',[]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
hold on; plot(xc(815)*1e-3,0,'or','markersize',15);
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'fontsize',15)
text(41,-220,'W (ms^-^1)','fontsize',15,'fontname','times');

%- shear square
ax4 = axes('position',[0.1 0.046 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
    [0:0.01:1,1:3:40]);
set(hh,'linewidth',0.1)
caxis([0 1]);%1
set(hh,'edgecolor','none');
set(gca,'xtick',[0:5:60]);
set(gca,'tickdir','out')
% cmocean('delta');
cmocean('-deep');
% colormap(ax4,'default');
% load MODIS_colorbar.mat;MODIS = colormap;clear colormap;
% colormap(ax4,(MODIS));
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
    [-0.2:0.1:0.2],'color','k');
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);

% set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
% set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
set(gca,'xticklabel',[{0} {0.5} {1.0} {1.5} {2.0} {2.5} {3.0} {3.5} {4.0} {4.5} {5.0} ]);
myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
z4 = colorbar('location','eastoutside');
set(z4,'position',[0.92 0.046 .02 0.22]); 
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
ylabel('Depth (m)','fontsize',15,'fontname','times')
xxx = xlabel('X (km)','fontsize',15,'fontname','times');
% get(xxx);
set(xxx,'position',get(xxx,'position')+[2.8 22 1.0]);
% set(xxx,'position',[0.0 .5 0.5]);
set(gca,'fontsize',15);
set(gcf,'color','w');
text(41,-220,'S^2 (10^-^3s^-^2)','fontsize',15,'fontname','times');

ax5 = axes('position',[0.66 0.33 .226 .164]);
U_interest = squeeze(U(815,1,1,:));
U_lowpass = lpass(U_interest,2/60,2,2);
% figure;
% plot((itv-itv(1))/3600,U_lowpass);
% hold on;plot((itv-itv(1))/3600,U_interest)
plot((itv-itv(1))/3600,U_lowpass,'k-','linewidth',2);
hold on; plot((itv(i)-itv(1))/3600,U_lowpass(i),'ro','markersize',15);
axis([0 25 0.2 2]);
ylabel('U(m/s)');xlabel('Time(hour)');grid on;
set(gca,'fontsize',10);
%------------------ time lag is too large for the boundary U
% t=linspace(0,24,length(itv));
% Period = 12.41667; %tidal period in hours
% ua=(-0.15*atan((-zc-200)/100)+0.5)*.8;% Along stream velocity
% for tt=1:length(t)
% Uvel(1,:,tt) = ua+.5*ua.*sin(2*pi/Period*t(tt));
% end
% plot(t,squeeze(Uvel(1,1,:)),'linewidth',2,'color','k');
% hold on; plot(t(i),squeeze(Uvel(1,1,i)),'ro');
% set(gcf,'color','w');
% xlabel('Time (hour)');ylabel('U (m/s)'); set(gca,'fontsize',15);grid on;
% ---------------------
drawnow 
 frame = getframe(gcf);
 if size(frame.cdata,1)~=926||size(frame.cdata,2)~=1065
%      figset(2,1);
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 31 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp)
set(gcf,'position',ps)
 frame = getframe(gcf);

 end
 
 writeVideo(temp,frame);
 close;
end
close(temp);