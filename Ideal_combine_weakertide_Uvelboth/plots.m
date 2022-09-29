%%% ITS = 39 flood tide; 33 ebb tide

clear all;load uvt2.mat;
addpath(genpath('/home/chihlun/MITgcm/mhlib'));
load('/home/chihlun/MITgcm/work/Bathymetry/topo.mat');
nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.15,122.418,nx);
y_new = linspace(24.29,24.9,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';

close all;
dim=size(U);
U_z=nan(dim(1),dim(2),dim(3),dim(4));
% x_range = 790:970;
x_range = 700:1100;

for i=[1:length(itv)]
for j=x_range %X-axis
ut=squeeze(U(j,:,:,i));
U_z(j,:,:,i)=mmderiv(-zc,ut')';
end
end

%- plot shear, V, W altogether 
%  U_z(U_z==0)=nan;T(T==0)=nan;U(U==0)=nan;W(W==0)=nan;V(V==0)=nan;
%  set(gca,'nextplot','replacechildren');
%  temp = VideoWriter('TUW_hd_time_wider1.avi');
%  temp.FrameRate = 2;
%  temp.Quality = 100;
%  open(temp);

% itv = ITS;
% for i=1:length(itv)
% for i=[1:382]
% for i=[1:length(itv)]

load storm_color.mat
%     for i=1:2:length(itv)
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
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[12.5:.05:28]);%caxi
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
% caxis([14 21]);
caxis([14 25]);
% z1 = colorbar('location','eastoutside');
% set(z1,'position',[0.92 0.76 .02 0.22]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-270 0]);
ylim([-330 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
set(gca,'xticklabel',[]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
% colormap(ax1,flipud(cbrewer('div','RdYlBu',100)));
colormap(ax1,(cm));
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
caxis([-.5 2]);
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
% ax4 = axes('position',[0.1 0.046 .8 .22]);
% [c,hh]=contourf(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
%     [0:0.01:1,1:3:40]);
% set(hh,'linewidth',0.1)
% caxis([0 1]);%1
% set(hh,'edgecolor','none');
% set(gca,'xtick',[0:5:60]);
% set(gca,'tickdir','out')
% % cmocean('delta');
% cmocean('-deep');
% % colormap(ax4,'default');
% % load MODIS_colorbar.mat;MODIS = colormap;clear colormap;
% % colormap(ax4,(MODIS));
% hold on;
% [c,hh]=contour(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
%     [-0.2:0.1:0.2],'color','k');
% set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
% 
% % set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end)))*1e-3);
% % set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
% set(gca,'xticklabel',[{0} {0.5} {1.0} {1.5} {2.0} {2.5} {3.0} {3.5} {4.0} {4.5} {5.0} ]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
% z4 = colorbar('location','eastoutside');
% set(z4,'position',[0.92 0.046 .02 0.22]); 
% hold on; 
% fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% % ylim([-270 0]);
% ylim([-330 0]);
% xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
% ylabel('Depth (m)','fontsize',15,'fontname','times')
% xxx = xlabel('X (km)','fontsize',15,'fontname','times');
% % get(xxx);
% set(xxx,'position',get(xxx,'position')+[2.8 22 1.0]);
% % set(xxx,'position',[0.0 .5 0.5]);
% 
% set(gca,'fontsize',15);
% set(gcf,'color','w');
% text(41,-220,'S^2 (10^-^3s^-^2)','fontsize',15,'fontname','times');

ax5 = axes('position',[0.66 0.33 .226 .164]);
U_interest = squeeze(U(815,1,1,:));
U_lowpass = lpass(U_interest,2/60,2,2);
% figure;
% plot((itv-itv(1))/3600,U_lowpass);
% hold on;plot((itv-itv(1))/3600,U_interest)
plot((itv-itv(1))/3600,U_lowpass,'k-','linewidth',2);
hold on; plot((itv(i)-itv(1))/3600,U_lowpass(i),'ro','markersize',15);
axis([0 25 1.5 2.5]);
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
%---------------------
% drawnow 
%  frame = getframe(gcf);
%  if size(frame.cdata,1)~=926||size(frame.cdata,2)~=1065
% %      figset(2,1);
% set(gcf,'units','centimeters','paperunits','centimeters')
% set(gcf,'PaperType','A4');
% pp=[0.63 0.9 31 28];
% ps=[0 0 pp(3)/1.1 pp(4)/1.1];
% set(gcf,'paperposition',pp)
% set(gcf,'position',ps)
%  frame = getframe(gcf);
% 
%  end
 
%  writeVideo(temp,frame);
%  close;
end
close(temp);



%% time series (same point as mooring G3: 122.3112, 24.6575)
%- xc(908): the data from this point
% [c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[16:.2:28]);%caxi
time = 0:2/60:2/60*(length(itv)-1);
time1 = 1; time2 = time1+length(itv)-1;

figure;figset(2,1);
% T(T==0)=nan;U(U==0)=nan;W(W==0)=nan;
% set(gcf,'units','centimeters','paperunits','centimeters')
% set(gcf,'PaperType','A4');
% pp=[0.63 0.9 27.5 28];
% ps=[0 0 pp(3)/1.1 pp(4)/1.1];
% set(gcf,'paperposition',pp)
% set(gcf,'position',ps)

ax1 = axes('position',[0.1 0.76 .8 .22]);
[c,hh]=contourf(time(time1:time2),zc,squeeze(T(908,1,:,time1:time2)),[16:.1:28]);
caxis([16.5 25]);
set(hh,'edgecolor','none');
% hold on;[c,hh]=contour(time(time1:time2),zc,squeeze(T(323,1,:,time1:time2)),[16:.2:22.5],'k');
hold on;[c,hh]=contour(time(time1:time2),zc,squeeze(T(908,1,:,time1:time2)),[16:.4:28],'k');

ylim([-200 -20]);
z1 = colorbar('location','eastoutside');
set(z1,'position',[0.92 0.76 .02 0.22]);
set(gca,'tickdir','out')
set(gca,'xtick',[time(time1):10:time(time2)]);
set(gca,'xticklabel',{[] [] [] [] [] [] [] [] [] [] [] [] []});
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
colormap(ax1,flipud(cbrewer('div','RdBu',100)));
% colormap(ax1,'default');


ax2 = axes('position',[0.1 0.522 .8 .22]);
[c,hh]=contourf(time(time1:time2),zc,squeeze(U(908,1,:,time1:time2)),[-1:.01:2]);
caxis([-.5 1.5]);
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(time(time1:time2),zc,squeeze(U(908,1,:,time1:time2)),[-1:.2:2],'k');
set(gca,'xtick',[time(time1):10:time(time2)]);
set(gca,'xticklabel',{[] [] [] [] [] [] [] [] [] [] [] [] []});
ylim([-200 -20]);
set(gca,'tickdir','out')
colormap(ax2,flipud(cbrewer('div','RdGy',100)));
z2 = colorbar('location','eastoutside');
set(z2,'position',[0.92 0.522 .02 0.22]);


ax3 = axes('position',[0.1 0.284 .8 .22]);
[c,hh]=contourf(time(time1:time2),zc,squeeze(W(908,1,:,time1:time2)),[-.5:.01:.5]);
caxis([-.4 .4]);
set(hh,'edgecolor','none');
z3 = colorbar('location','eastoutside');
set(z3,'position',[0.92 0.284 .02 0.22]);
set(gca,'tickdir','out')
ylim([-200 -20]);
colormap(ax3,flipud(cbrewer('div','RdGy',100)));
% hold on;
% [c,hh]=contour(time(time1:time2),zc,squeeze(W(825,1,:,time1:time2)),'k');

% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((W((x_range),1,:,i)))',[-1:0.5:2],'k');
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
set(gca,'xtick',[time(time1):1:time(time2)]);
set(gca,'xticklabel',{'0','30','60','90','120','150','180','210','240','270','300','330','360'});
xlabel('minute','fontsize',10);
ylabel('D');
set(gcf,'color','w');
print -dpdf -r600 Mooring_time_series_comparison.pdf



%% streamline 

stream = VideoWriter('Streamline_T.avi');
stream.FrameRate = 3;
stream.Quality = 100;
open(stream);
nx = 1377;
[LON,LAT] = meshgrid(x,y);
p1 = [122.11,24.2]; p2 = [122.365,24.775]; %2 points where the domain is
x_new = linspace(p1(1),p2(1),nx);
y_new = linspace(p1(2),p2(2),nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';
%%% find the mooring spot
% R = 6371*1e3;
% a = sind((p2(2)-p1(2))/2).^2 + cosd(p1(2))*cosd(p2(2)*sind((p2(1)-p1(1))/2)^2);
% c = 2*atan2d(sqrt(a),sqrt(1-a));
% d = R*c;


x_range = 790:970;
for i=1:10
figure;figset(2,1);
ax1 = axes('position',[0.1 0.55 .8 .4]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[14:.1:28.5]);%caxi
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
set(z1,'position',[0.905 0.55 .015 0.4]);
ylabel(z1,'Temperature(^oC)');
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
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
tit = title(sprintf('Time=%.1f min',2*(i-1)));
set(tit,'fontsize',10)
text(41.75,-220,'T (^oC)','fontsize',15,'fontname','times');


%-streamline
% figure;figset(2,1);
xrange = [690]; zrange = 1:112;
[xx,zz] = meshgrid(xc,zc);
[sx,sz] = meshgrid(xc(xrange),zc(zrange));
ax2 = axes('position',[0.1 0.1 .8 .4]);
XZ1 = streamline(xx,zz,squeeze(U(:,1,:,i))',squeeze(W(:,1,:,i))',sx,sz,[0.1,20000]);
set(XZ1,'color','k');


hold on; plot(xc(1099),zc(1),'r^');
hold on; plot(ones(length(zc),1)*xc(1099),zc(:),'r');
hold on;
fill(xc,zb,[190 190 190]/225); % topo black
set(gca,'xtick',(xc(x_range(1)):500:xc(x_range(end))));
set(gca,'xticklabel',{'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'});
ylim([-250 0])
xlim([xc(x_range(1)) xc(x_range(end))]);
set(gca,'tickdir','out')


xlabel('km','fontsize',10);ylabel('D(m)','fontsize',10);
set(gcf,'color','w');
% drawnow 
%  frame = getframe(gcf);
%  if size(frame.cdata,1)~=955||size(frame.cdata,2)~=673
%      figset(2,1);
%  frame = getframe(gcf);
%  end
% 
% writeVideo(stream,frame);
% close;

end
close(stream);

%% 
close all;

dim=size(U);
U_z=nan(dim(1),dim(2),dim(3),dim(4));
x_range = 21*i+720:3:21*i+742;
z_range = 1:2:110;

for i=315
for j=1:dim(1)  %X-axis
ut=squeeze(U(j,:,:,i));
U_z(j,:,:,i)=mmderiv(-zc,ut')';
end
end
i=4;
figure;figset(2,1)

[c,hh]=contourf(xc(x_range),zc(z_range),1e3*(squeeze(U_z(x_range,1,z_range,315)))',[-160:20]);
% set(hh,'linewidth',0.1)
caxis([-80 20]);%1
set(hh,'edgecolor','none');
set(gca,'tickdir','out')
colormap(flipud(cbrewer('seq','PuRd',100)));
cmocean('-haline') ;
hold on;
[c,hh]=contour(xc(x_range),zc(z_range),1e3*(squeeze(U_z(x_range,1,z_range,315)))',[0 0],'r');
set(hh,'linewidth',2)

% colormap(flipud(brewermap([],'Spectral')));

% colormap default
hold on;

% z_range = 10:112;
hhh=quiver(xc(x_range),zc(z_range),28*squeeze(U(x_range,1,z_range,315))',...
    28*squeeze(W(x_range,1,z_range,315))',...
    'color','k');
axis equal;
set(hhh,'autoscale','off');
set(hhh,'maxheadsize',0.1,'linewidth',1);
set(gcf,'color','w');
set(gca,'xtick',[xc(x_range(1)):50:xc(x_range(end))]);
set(gca,'xticklabel',{'0','50','100','150','200','250','300','350','400'});
hold on;
fill(xc,zb,[190 190 190]/225); % topo black
% ylim([-330 0])
xlim([xc(x_range(1)) xc(x_range(end))]);

aaa = colorbar;ylabel(aaa,'\partialU/\partialz (s^-^1)','fontsize',15);
xlabel('X(m)');ylabel('D(m)');
set(gca,'fontsize',20);
ylim([-300,0]);
hhh=quiver(xc(x_range(2)),-190,25*1,0,'color','k');
set(hhh,'autoscale','off');
set(hhh,'maxheadsize',0.8,'linewidth',1);
text(xc(x_range(2)),-200,'1ms^-^1','fontsize',12);




%% Fr number estimation Fr = (zou0*U^2)/(g*D^2*zou0_z)
load('/home/chihlun/Data_base/200mData/topo.mat');
nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.11,122.365,nx);
y_new = linspace(24.2,24.775,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';
clear dens_z
% potential density
dim=size(U);
time = 250:320;
pdens=nan(dim(1),dim(2),dim(3),length(time));

S=35; % constant salinity
% potential density
%for i=1:time
for i=time
for j=1:dim(3)
pdens(:,:,j,i)=sw_pden(S*ones(dim(1),1),T(:,1,j,i),-zc(j),0);
end
end
for i=time
for j=1:dim(1)
dens_z(j,1,:,i)=mmderiv(-zc,squeeze(pdens(j,1,:,i))')';
end
end
dens_zt = mean(dens_z(:,1,:,time),4);
%- reorder 
% pdens_reorder = sort(squeeze(pdens(:,1,:,time)),2);
% pdens_purt = squeeze(pdens(:,1,:,time))-pdens_reorder;
% % for i=time
% for j=1:dim(1)
% dens_z(j,:)=mmderiv(-zc,squeeze(pdens_purt(j,:))')';
% end
% end
%%%%%%%%%%%%%%%%%%%%  plotting the vertical density gradient to compare
%%%%%%%%%%%%%%%%%%%%  with echo intensity
% x_range = 800:860;
% % itv = ITS;
% % for i=1:length(itv)
% figure;
% [c,hh]=contourf(xc(x_range)*1e-3,zc(1:end),squeeze(dens_z((x_range),1:end))',...
%     [-.1:0.01:0,0.00005:0.00005:0.1,0.2:0.2:1]);%caxi
% 
% % [c,hh]=contourf(xc(x_range)*1e-3,zc(18:end),squeeze(dens_z((x_range),1,18:end,i))',...
% %     [-0.01:.0001:0.01]);%caxi
% colormap((flipud(cbrewer('seq','YlGnBu',100))));
% caxis([0 0.02]);
% % set(gca,'color','k');
% % hold on;
% % [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[22:.4:28],'color',[40 40 40]/225 );
% % hold on;
% % [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[14:.4:16.9],'color',[40 40 40]/225);
% set(hh,'edgecolor','none');
% 
% 
% z1 = colorbar('location','eastoutside');
% set(gca,'tickdir','out')
% hold on; 
% fill(xc*1e-3,zb,[190 190 190]/225); % topo black
% ylim([-250 0])
% xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
% set(gca,'xtick',linspace(xc(x_range(1)),xc(x_range(end)),4)*1e-3);
% set(gca,'xticklabel',{'0','1.1','2.2','3.3'});
% 
% set(gca,'fontsize',15);
% % text(42.5,-200,'\partial','fontsize',15,'fontname','times');
% set(gcf,'color','w')


%%%%%%%%%%%%%%%%%%%%
%- estimation of the Froude number 
% dens_zt = (dens_z(:,1,:,1));
dens_z_avg_tz = mean(dens_zt(:,1:18),2); %average the upper 50 m 
zou0 = 1026;
g = 9.81;
D = zb;
D(D<zc(end))=zc(end);

UU=[];
for i = time
UU = [UU,mean(squeeze(U(:,1,1:18,i)),2)]; %average the upper 50 m 
end
UU = mean(UU,2);

% UU = mean(UU,2);
Fr = sqrt((zou0*UU.^2)./(g*D.^2.*dens_z_avg_tz));
figure;plot(xc(790:850)*1e-3,Fr(790:850),'-','linewidth',2,'color','b');
% hold on; plot(linspace(xc(700),xc(970),100),1*ones(1,100),'-r-');
xlabel('X(km)','fontsize',20);ylabel('Fr','fontsize',20);
set(gcf,'color','w');
% xlim([41.5 43]);
set(gca,'fontsize',15,'tickdir','out');
set(gca,'xtick',[xc(790):500:xc(850)]*1e-3);
myticks('x','y',1e-3*[xc(790):100:xc(850)],.98);
ylim([0.25 0.285]);
set(gca,'xticklabel',{'0','0.5','1.0'});
set(gcf,'color','w');
grid on
Dc = zb(find(Fr>1)); % the critical depth 
% my own definition of the upstream mean flow
U_upstream = mean(UU(300:700)); 
x_range = 1000:1100;
EE = -zb(x_range)+UU(x_range).^2/(2*9.81);
figure;plot(EE,linspace(0,max(-zb(x_range)),length(x_range)));


%% model configuration 

load UVWT.mat ;
addpath(genpath('/home/chihlun/MITgcm_c65r/mhlib'));
load('/home/chihlun/Data_base/200mData/topo.mat');
nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.15,122.418,nx);
y_new = linspace(24.29,24.9,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';

x_range = 1:nx;
figure;
%%%%%%%%%%%% figure size
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 27.5 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
%%%%%%%%%%%%

ax1 = axes('position',[0.1 0.58 .8 .4]);
[c,hh]=contourf(xc*1e-3,zc,squeeze(T(:,1,:,1))',[13:.1:28]);
% caxis([16.5 25]);
set(hh,'edgecolor','none');
z1 = colorbar('location','eastoutside');
set(z1,'position',[0.92 0.58 .02 0.4]);
set(gca,'tickdir','out');
% set(gca,'xtick',[time(time1):30:time(time2)]);
% set(gca,'xticklabel',{[] [] [] [] [] [] [] [] [] [] [] [] []});
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
colormap(ax1,flipud(cbrewer('div','RdYlBu',100)));
hold on;
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-336 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
% hold on;
% plot(ones(1,112)*xc(614)*1e-3,zc(1:112),'color','r','linewidth',2);
% plot(ones(1,112)*xc(845)*1e-3,zc(1:112),'color','r','linewidth',2);

set(gca,'xtick',[]);set(gca,'ytick',[]);
% xlabel('X = 72.89 km','fontsize',15,'fontweight','bold');
% ylabel('Z = 334.5 m','fontsize',15,'fontweight','bold');
% set(gca,'ytick',[-336,0]);
% set(gca,'xtick',linspace(xc(x_range(1)),xc(x_range(end)),11)*1e-3);
% set(gca,'xticklabel',{[],[],[],[],[],[],[],[],[],[]});
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
% colormap(ax1,flipud(cbrewer('div','RdBu',100)));
% colormap(ax1,'default');
% axis equal 
% text(60,-250,'\Deltax = 50 m','fontsize',15,'fontweight','bold','color','w');
% text(60,-270,'\Deltaz = 3 m','fontsize',15,'fontweight','bold','color','w');
% text(60,-290,'\Deltat = 3 s','fontsize',15,'fontweight','bold','color','w');
hold on; 
xx=quiver(.1,zc(1:5:end),30*squeeze(U(1,1,1:5:end,1)),0*ones(23,1),'color','k','linewidth',.3);
xx.MaxHeadSize=0.3;
set(xx,'autoscale','off');set(gca,'fontsize',15)

% ylim([-336 5])
ax2 = axes('position',[0.048 0.1 .18 .4]);
Temper = plot(squeeze(T(1,1,:,1)),zc,'-ok');
set(Temper,'linewidth',1.2,'markersize',3);
% xlabel('T (^oC)','fontsize',15,'fontweight','bold');
% ylabel('Z (m)','fontsize',15,'fontweight','bold');
set(gca,'tickdir','out')
set(gca,'fontsize',15)
grid minor

S = 35;
ax3 = axes('position',[0.263 0.1 .18 .4]);
Sal = plot(squeeze(S*ones(length(zc),1)),zc,'o-r');
set(Sal,'linewidth',1.2,'markersize',3);
xlabel('S = 35 psu','fontsize',15,'fontweight','bold');
set(gca,'yticklabel',[],'tickdir','out')
set(gca,'fontsize',15)

grid minor

ax4 = axes('position',[0.483 0.1 .18 .4]);
UUU = plot(squeeze(U(1,1,:,1)),zc,'o-b');
set(UUU,'linewidth',1.2,'markersize',3);
xlim([0 .8])
set(gca,'xtick',[0:.2:1])
% xlabel('U (ms^-^1)','fontsize',15,'fontweight','bold');
set(gca,'yticklabel',[],'tickdir','out')
set(gca,'fontsize',15)

grid minor

ax5 = axes('position',[0.77 0.1 .22 .4]);
dx = [340:-2:20,20*ones(1,1075),20:2:300]; %dy=dx; 
xf(1)=0;
for i=1:nx
xf(i+1)=xf(i)+dx(i);
end
xc=(xf(1:nx)+xf(2:nx+1))/2;

UUU = plot(xc*1e-3,dx,'o-k');
set(UUU,'linewidth',1.2,'markersize',3);
hold on;
plot(ones(1,351)*xc(614)*1e-3,0:350,'color','r','linewidth',2);
plot(ones(1,351)*xc(845)*1e-3,0:350,'color','r','linewidth',2);
set(gca,'xtick',[0:15:60,72]);
set(gca,'ytick',[0,20,50:50:350]);
% ylabel('Horizontal Resolution \Deltax (m)')
% xlabel('X (km)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15)
grid minor


set(gcf,'color','w');

print('-fillpage','Model_setup','-painters','-r1000','-dpdf');

%% what kind of phenomenon is this? 

% zc(18)=52.5, zc(43)=127.5, xc(825)=42250, xc(888)=43510, xc(875)=43250
% ezfft U
close all
clear de_U 
%figure; plot(Lon,Lat);
init = 1;
final = 324;
% final = final(end);
% de_U=squeeze(T(1099,1,48,init:final));
de_U=detrend(squeeze(T(875,1,43,init:final)));
%figure;plot(de_U);

figure;figset(2,1);
subplot(2,1,1);
% plot(ITS(init-(init-1):final-(init-1))/1200,de_U);
[W,E]=ezfft(itv(init-(init-1):final-(init-1))/60,de_U,'disp'); %60 min=1 hour
Frequency=W/2/pi;
% set(gca,'fontsize',14);
subplot(2,1,2);
plot(1./(W./(2*pi)),E*1e3,'linewidth',1);
xlim([0 60]);
myticks('x','y',[0:2:60],.8);
set(gcf,'color','w');
% set(gca,'fontsize',14);

%set(gca,'xtick',[0:2:60]);
%set(gca,'xticklabel',{'0',[],[],[],[],'10',[],[],[],[],'20',...
    %[],[],[],[],'30',[],[],[],[],'40',[],[],[],[],'50',[],[],[],[],'60'});

xlabel('Period (Minute)','fontsize',15,'fontweight','bold');
ylabel('Energy Density (10^-^3 m^2/Hours)','fontsize',15,'fontweight','bold');
print('-dpdf','Energy Density Spectrum U.pdf');







%% - 1) Internal Solitary Waves (ISW)
clear all;load uvt2.mat;
addpath(genpath('/home/chihlun/MITgcm_c65r/mhlib'));
load('/home/chihlun/Data_base/200mData/topo.mat');
nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.15,122.418,nx);
y_new = linspace(24.29,24.9,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';


x_range = 790:900;
i=10;
figure;
%%%%%%%%%%%% figure size
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 27.5 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp)
set(gcf,'position',ps)
%%%%%%%%%%%%

%- T
ax1 = axes('position',[0.08 0.55 .4 .4]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[14:.1:28.5]);%caxi
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[15:.3:28],'color',[40 40 40]/225);
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[22 22],'color','r');
set(hh,'linewidth',2);
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[17 17],'color','g');
set(hh,'linewidth',2);
caxis([15 28]);

z1 = colorbar('location','southoutside');
set(z1,'position',[0.08 0.48 .4 0.02]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',linspace(xc(x_range(1)),xc(x_range(end)),5)*1e-3);
set(gca,'xticklabel',{'0','0.55','1.1','1.65','2.2'});
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
colormap(ax1,flipud(cbrewer('div','RdYlBu',100)));
% colormap(ax1,'default');
% colormap(ax1,flipud(brewermap([],'Spectral')))
% load MODIS_colorbar.mat;MODIS = colormap;clear colormap;
% colormap(ax1,(MODIS));
% xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'fontsize',15);
% tit = title(sprintf('Time=%.1f min',2*(i-1)));
% set(tit,'fontsize',10)
text(41.75,-220,'T (^oC)','fontsize',15,'fontname','times');

%- W
ax2 = axes('position',[0.53 0.55 .4 .4]);
% ax3 = axes('position',[0.1 0.284 .8 .22]);
 [c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((W(x_range,1,:,i)))',[-.5:0.005:.5]);
 caxis([-.4 .4]);
 set(hh,'edgecolor','none');
 z2 = colorbar('location','southoutside');
set(z2,'position',[0.53 0.48 .4 0.02]);
set(gca,'tickdir','out')
colormap(ax2,flipud(cbrewer('div','RdGy',100)));
colormap(ax2,'jet');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,i))',[22 22],'color','r');
set(hh,'linewidth',2);

% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((W((x_range),1,:,i)))',[-1:0.5:2],'k');
set(gca,'xtick',linspace(xc(x_range(1)),xc(x_range(end)),5)*1e-3);
set(gca,'xticklabel',{'0','0.55','1.1','1.65','2.2'});
set(gca,'yticklabel',[]);
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'fontsize',15)
text(41.75,-220,'W (ms^-^1)','fontsize',15,'fontname','times');
set(gcf,'color','w')
print('-fillpage','ISW','-painters','-r1000','-dpdf');

%- estimate the ISW's vertical displacement and propagating speed 

% C = (du/dt) / (dw/dz)


dim=size(W);
W_z=nan(dim(1),dim(2),dim(3),dim(4));
U_t=nan(dim(1),dim(2),dim(3),151);
pos = 850;

for i=1:dim(4)
for j= pos %this is where the ISW occurs
wt=squeeze(W(j,:,:,i));
W_z(j,:,:,i)=mmderiv(-zc,wt')';
end
end
W_z=squeeze(W_z(pos,1,:,:));
time = [0:150]*120;

for i=1:dim(3)
    for j=pos
        ut=squeeze(U(j,:,i,1:151));
        U_t(j,:,i,:)=mmderiv(time,ut');
    end
end
U_t=squeeze(U_t(pos,1,:,:));

figure;figset(1,2);
for i=1:150
    scatter(W_z(1:50,i)*1e3,U_t(1:50,i)*1e3,10,zc(1:50),'filled');
    hold on;
end
colorbar;
caxis([-150 0])

hold on;
p = polyfit(W_z(1:50,1:150)*1e3,U_t(1:50,1:150)*1e3,1);
f = polyval(p,W_z(1:50,1:150)*1e3); 
colormap(flipud(brewermap([],'YlGnBu')));
set(gca,'color','w');
prop = plot(W_z(1:50,1:150)*1e3,f) ;
set(prop,'color','r','linestyle','-','linewidth',2);

axis image;
box on;
set(gcf,'color','w');
axis([-12 8 -6 8]);
set(gca,'fontsize',15,'tickdir','out');
myticks('x','y',-15:1:10,.8);
myticks('y','y',-5:1:7,.8);
xlabel('\partialw/\partialz (10^-^3 s^-^1)')
ylabel('\partialu/\partialt (10^-^3 ms^-^2)')


%% 2) hydraulic jump and Shear Instability
% 46 51 58 

i = [46,51,58];
x_range = 790:970;

figure;
%%%%%%%%%%%% figure size
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 27.5 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp)
set(gcf,'position',ps)
%%%%%%%%%%%%

%- T
ax1 = axes('position',[0.1 0.74 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,46))',[14:.1:28.5]);%caxi
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,46))',[15:.3:28],'color',[40 40 40]/225);
caxis([15 27]);
% z1 = colorbar('location','eastoutside');
% set(z1,'position',[0.92 0.74 .02 0.22]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',[xc(x_range(1)):500:xc(x_range(end))]*1e-3);
myticks('x','y',[xc(x_range(1)):100:xc(x_range(end))]*1e-3,.8,'popupid');
myticks('y','y',-300:10:0,.8,'popupid');
set(gca,'xticklabel',[]);
% set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
colormap(ax1,flipud(cbrewer('div','RdYlBu',100)));
set(gca,'fontsize',15);
% tit = title(sprintf('Time = t_0 (min)'));
% set(tit,'fontsize',10)
text(41.63,-220,'t_0 (min)','fontsize',15,'fontname','times');

ax2 = axes('position',[0.1 0.48 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,51))',[14:.1:28.5]);%caxi
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,51))',[15:.3:28],'color',[40 40 40]/225);
caxis([15 27]);
% z1 = colorbar('location','eastoutside');
% set(z1,'position',[0.92 0.76 .02 0.22]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',[xc(x_range(1)):500:xc(x_range(end))]*1e-3);
myticks('x','y',[xc(x_range(1)):100:xc(x_range(end))]*1e-3,.8,'popupid');
myticks('y','y',-300:10:0,.8,'popupid');
set(gca,'xticklabel',[]);
% set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
colormap(ax2,flipud(cbrewer('div','RdYlBu',100)));
set(gca,'fontsize',15);
% tit = title(sprintf('Time = t_0+10 (min)'));
% set(tit,'fontsize',10)
text(41.63,-220,'t_0+10 (min)','fontsize',15,'fontname','times');


ax3 = axes('position',[0.1 0.22 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,58))',[14:.1:28.5]);%caxi
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,58))',[15:.3:28],'color',[40 40 40]/225);
caxis([15 27]);
z3 = colorbar('location','eastoutside');
set(z3,'position',[0.92 0.22 .02 0.22]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',[xc(x_range(1)):500:xc(x_range(end))]*1e-3);
myticks('x','y',[xc(x_range(1)):100:xc(x_range(end))]*1e-3,.8,'popupid');
myticks('y','y',-300:10:0,.8,'popupid');
% set(gca,'xticklabel',[]);
set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
colormap(ax3,flipud(cbrewer('div','RdYlBu',100)));
set(gca,'fontsize',15);
% tit = title(sprintf('Time = t_0+24 (min)'));
% set(tit,'fontsize',10)
text(41.63,-220,'t_0+24 (min)','fontsize',15,'fontname','times');
% text(41.63,-220,'T(^oC)','fontsize',15,'fontname','times');


xlabel('x(km)','fontsize',15);ylabel('D(m)','fontsize',15);
set(gcf,'color','w');
% print('KH','-painters','-r1000','-dpdf');
print('-fillpage','KH','-painters','-r600','-dpdf');


%% Richardson Number and dissipation rate

load('/home/chihlun/Data_base/200mData/topo.mat');

nx = 1377;
[LON,LAT] = meshgrid(x,y);
x_new = linspace(122.15,122.418,nx);
y_new = linspace(24.29,24.9,nx);
zb = interp2(LON,LAT,z,x_new,y_new,'linear')';
% creating the vertical shear S
dim=size(U);
U_z=nan(dim(1),dim(2),dim(3),dim(4));

for i=58
for j=1:dim(1)  %X-axis
ut=squeeze(U(j,:,:,i));
U_z(j,:,:,i)=mmderiv(-zc,ut')';
end
end
S2 = U_z.^2;

%creating the buoyancy frequency N2
clear dens_z
% potential density
dim=size(U);
time = 58;
pdens=nan(dim(1),dim(2),dim(3),length(time));

S=35; % constant salinity
% potential density
%for i=1:time
for i=time
for j=1:dim(3)
pdens(:,:,j,i)=sw_pden(S*ones(dim(1),1),T(:,1,j,i),-zc(j),0);
end
end
% for i=time
% for j=1:dim(1)
% dens_z(j,1,:,i)=mmderiv(-zc,squeeze(pdens(j,1,:,i))')';
% end
% end
% dens_zt = mean(dens_z(:,1,:,time),4);

% reorder the density
pdens_reorder = sort(squeeze(pdens(:,1,:,time)),2);
pdens_purt = squeeze(pdens(:,1,:,time))-pdens_reorder;
for i=time
for j=1:dim(1)
dens_z(j,:)=mmderiv(-zc,squeeze(pdens_reorder(j,:))')';
% dens_z(j,:)=mmderiv(-zc,squeeze(pdens(j,1,:,i))')';
end
end
% some parameters
g = 9.81; zou0 = 1026;
N2 = g/zou0*dens_z;

Ri = N2./squeeze(S2(:,1,:,58)); % richardson number 
% Ri(Ri>100|Ri<-100)=nan;
x_range = 790:970;


figure;
%%%%%%%%%%%% figure size
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
pp=[0.63 0.9 27.5 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp)
set(gcf,'position',ps)
%%%%%%%%%%%%

ax1 = axes('position',[0.1 0.76 .8 .22]);

% hh=pcolor(xc(x_range)*1e-3,zc(1:end),1./Ri(x_range,:)');
[cc,hh] = contourf(xc(x_range)*1e-3,zc(1:end),1./Ri(x_range,:)',[0:0.1:50]);
set(hh,'edgecolor','none');
caxis([4 25]);
% colormap(ax1,flipud(cbrewer('div','RdYlBu',100)));
colormap(ax1,'hot');

z1 = colorbar('location','eastoutside');
set(z1,'position',[0.92 0.76 .02 0.22]);
set(z1,'ytick',[4,10 15 20 25]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',[xc(x_range(1)):500:xc(x_range(end))]*1e-3);
% set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
set(gca,'xticklabel',[]);
% 
set(gca,'fontsize',15);
text(41.75,-200,'Ri^-^1','fontsize',15,'fontname','times');
set(gcf,'color','w');
myticks('x','y',[xc(x_range(1)):100:xc(x_range(end))]*1e-3,.8,'popupid');
myticks('y','y',-300:10:0,.8,'popupid');

%%%%%%%%%%%%%%%%%%%%%%%% dissipation rate

%- First part is to calculate the ensemble average
init = 1;
final = 324;
ensavg_U=mean(U(:,:,:,init:final),4);
ensavg_W=mean(W(:,:,:,init:final),4);
% ensavg_T=mean(T(:,:,:,init:final),4);
for i=init:final
U_prime(:,:,:,i-(init-1))=U(:,:,:,i-(init-1))-ensavg_U;
W_prime(:,:,:,i-(init-1))=W(:,:,:,i-(init-1))-ensavg_W;
end

% calculate dissipation rate epsilon
for tt=58
for i=1:length(zc)
du_dx(:,i,tt) = mmderiv(xc,squeeze(U_prime(:,:,i,tt)));
dw_dx(:,i,tt) = mmderiv(xc,squeeze(W_prime(:,:,i,tt)));
end
end

for tt=58
for i=1:length(xc)
du_dz(i,:,tt) = mmderiv(zc',squeeze(U_prime(i,:,:,tt)));
end
end
nu = 3;
epsilon = nu*.5.*((du_dx+du_dx).^2+(dw_dx+du_dz).^2);

x_range = 790:970;
ax2 = axes('position',[0.1 0.522 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,1e3*epsilon(x_range,:,tt)',[0:0.001:.5,.6:.1:1,2:30]);
caxis([0 .5]);
set(hh,'edgecolor','none');
z2 = colorbar('location','eastoutside');
set(z2,'position',[0.92 0.522 .02 0.22]);
set(gca,'tickdir','out')
% colormap(ax2,flipud(cbrewer('div','RdGy',100)));
colormap(ax2,'hot');
% set(gca,'xticklabel',{[],[],[],[],[],[],[],[],[],[]});
% myticks('x','x',linspace(xc(x_range(1)),xc(x_range(end)),10)*1e-3,.8,'popupid');
hold on; 
fill(xc*1e-3,zb,[190 190 190]/225); % topo black
ylim([-270 0])
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',[xc(x_range(1)):500:xc(x_range(end))]*1e-3);
set(gca,'xticklabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
set(gca,'yticklabel',{'-250','-200','-150','-100','-50',[]});
set(gcf,'color','w');
% text(41.75,-200,'Ri^-^1','fontsize',15,'fontname','times');
text(41.6,-200,'\epsilon(10^-^3xm^2s^-^3)','fontsize',15,'fontname','times');
set(gca,'fontsize',15);
myticks('x','y',[xc(x_range(1)):100:xc(x_range(end))]*1e-3,.8,'popupid');
myticks('y','y',-300:10:0,.8,'popupid');

xlabel('X(km)','fontsize',15);ylabel('D(m)','fontsize',15);
print('-fillpage','Ri_epsilon','-painters','-r600','-dpdf');














%%%%%%%%%%%%%%%%%%%%
%- estimation of the Froude number 
% dens_zt = (dens_z(:,1,:,1));
dens_z_avg_tz = mean(dens_zt(:,1:18),2); %average the upper 50 m 
zou0 = 1026;
g = 9.81;
D = zb;
D(D<zc(end))=zc(end);

UU=[];
for i = time
UU = [UU,mean(squeeze(U(:,1,1:18,i)),2)]; %average the upper 50 m 
end
UU = mean(UU,2);

% UU = mean(UU,2);
Fr = sqrt((zou0*UU.^2)./(g*D.^2.*dens_z_avg_tz));
figure;plot(xc(790:850)*1e-3,Fr(790:850),'-','linewidth',2,'color','b');
% hold on; plot(linspace(xc(700),xc(970),100),1*ones(1,100),'-r-');
xlabel('X(km)','fontsize',20);ylabel('Fr','fontsize',20);
set(gcf,'color','w');
% xlim([41.5 43]);
set(gca,'fontsize',15,'tickdir','out');
set(gca,'xtick',[xc(790):500:xc(850)]*1e-3);
myticks('x','y',1e-3*[xc(790):100:xc(850)],.98);
ylim([0.25 0.285]);
set(gca,'xticklabel',{'0','0.5','1.0'});
set(gcf,'color','w');
grid on
Dc = zb(find(Fr>1)); % the critical depth 
% my own definition of the upstream mean flow
U_upstream = mean(UU(300:700)); 
x_range = 1000:1100;
EE = -zb(x_range)+UU(x_range).^2/(2*9.81);
figure;plot(EE,linspace(0,max(-zb(x_range)),length(x_range)));















