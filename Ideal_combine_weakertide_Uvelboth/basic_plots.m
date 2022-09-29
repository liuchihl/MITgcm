addpath(genpath('D:/NTU2015-2019/NTU resource/TOOLS'));
clear all;%close all;%load uvt2_tide.mat;
load('uvt_all.mat');
load storm_color.mat % the colorbar for temperature contour
cm_U = load('NCV_rainbow2.rgb');
cm_Uz = load('MPL_gnuplot.rgb');


xc = data.X; yc = data.Y; zc = data.Z; itv = data.T;
U = data.U(:,length(yc)/2,:,:); V = data.V(:,length(yc)/2,:,:);
W = data.W(:,length(yc)/2,:,:); T = data.Temp(:,length(yc)/2,:,:);
time = 1:length(itv);
U_interest = squeeze(U(771,1,1,time));
U_lowpass = lpass(U_interest,2/60,2,2);
xticks_new = [22637:1000:xc(771):1000:42637]*1e-3;
x_range = 600:970;
clear data
%-1 shear
dim=size(U);
U_z=nan(dim(1),dim(2),dim(3));
for i=time
for j=x_range %X-axis
ut=squeeze(U(j,:,:,i));
U_z(j,:,:,i-time(1)+1)=mmderiv(zc,ut')';
end
end

%%
temp = VideoWriter('Ilan_3D.avi');
temp.FrameRate = 7;
temp.Quality = 100;
open(temp);


% close all;
h = figure;
%%%%%%%%%%%% figure size
set(gcf,'units','centimeters','paperunits','centimeters')
set(gcf,'PaperType','A4');
% pp=[0.63 0.9 27.5 28];
pp=[0.63 0.9 25 28];
ps=[0 0 pp(3)/1.1 pp(4)/1.1];
set(gcf,'paperposition',pp);
set(gcf,'position',ps);
%%%%%%%%%%%%

x_range = 700:970;
% x_range = 850:970;
% x_range = 1:1377;
for i=1:length(itv)
%- T
ax1 = axes('position',[0.1 0.76 .8 .22]);
% ax1 = axes('position',[0.1 0.45 .8 .5]);
          
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,time(i)))',linspace(10.8,28,262));%caxi
% figure;pcolor(xc*1e-3,zc,(squeeze(T(,1,:,i))'));
% shading flat
set(hh,'edgecolor','none');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze(T((x_range),1,:,time(i)))',[12:0.3:28],'color',[45 45 44]/225);
set(hh,'linewidth',.5);
caxis([12 26]);

z1 = colorbar('location','eastoutside');
set(z1,'position',[0.92 0.76 .02 0.22]);
% set(z1,'position',[0.92 0.45 .02 0.5]);
set(gca,'tickdir','out')
hold on; 
fill(xc*1e-3,zb3D,[190 190 190]/225); % topo black
ylim([-300 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',xticks_new);
set(gca,'ytick',-300:100:0);
set(gca,'xticklabel',[]);
% myticks('x','y',xticks_new(1):.2:xticks_new(end),.8,'popupid');
% myticks('y','y',-300:20:0,.8,'popupid');
colormap(ax1,flipud(cm));
set(gca,'fontsize',15);
% tit = title(sprintf('Time=%.1f minutes',2*(time(i)-time(1))));
tit = title(sprintf('$$Time = t_o+%.0f (min)$$',2*(time(i)-time(1))),'interpreter','latex');

set(tit,'fontsize',10,'fontweight','bold');
text(41,-220,'T (^oC)','fontsize',18,'fontname','arial');
hold off;

%- U
ax2 = axes('position',[0.1 0.522 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((U(x_range,1,:,time(i))))',[-1:0.01:2]);
% caxis([-.5 2.0]);
caxis([-.5 2.0]);
set(hh,'edgecolor','none');
hold on;
% [c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((U((x_range),1,:,i)))',[-1:0.2:2],'k');
[c,hh]=contour(xc(x_range)*1e-3,zc,squeeze((U((x_range),1,:,time(i))))',[0 0],'r');
z2 = colorbar('location','eastoutside');
set(z2,'position',[0.92 0.522 .02 0.22]);
colormap(ax2,cm_U/255);
hold on;fill(xc*1e-3,zb3D,[190 190 190]/225); % topo black
annotation('line',[.92 .939],[.567 .567],'color','r','linewidth',2);
ylim([-300 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',xticks_new);
set(gca,'ytick',-300:100:0);
set(gca,'xticklabel',[]);set(gca,'fontsize',15);
myticks('x','y',xticks_new(1):.2:xticks_new(end),.8,'popupid');
myticks('y','y',-300:20:0,.8,'popupid');
set(gca,'tickdir','out')
set(gca,'fontsize',15)
text(41,-220,'U (ms^-^1)','fontsize',18,'fontname','arial');

%- W
ax3 = axes('position',[0.1 0.284 .8 .22]);
 [c,hh]=contourf(xc(x_range)*1e-3,zc,squeeze((W(x_range,1,:,time(i))))',[-1:0.005:1]);
 caxis([-.5 .5]);
 set(hh,'edgecolor','none');
z3 = colorbar('location','eastoutside');
set(z3,'position',[0.92 0.284 .02 0.22]);
colormap(ax3,cm_U/255);
hold on;fill(xc*1e-3,zb3D,[190 190 190]/225); % topo black
annotation('line',[.92 .939],[.567 .567],'color','r','linewidth',2);
ylim([-300 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
set(gca,'xtick',xticks_new);
set(gca,'ytick',-300:100:0);
set(gca,'xticklabel',[]);set(gca,'fontsize',15);
myticks('x','y',xticks_new(1):.2:xticks_new(end),.8,'popupid');
myticks('y','y',-300:20:0,.8,'popupid');
set(gca,'tickdir','out')
set(gca,'fontsize',15)
text(41,-220,'W (ms^-^1)','fontsize',18,'fontname','arial');

%- shear square
ax4 = axes('position',[0.1 0.046 .8 .22]);
[c,hh]=contourf(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
    [0:0.01:2,1:2:40]);
set(hh,'linewidth',0.1)
caxis([0 1]);%1
set(hh,'edgecolor','none');
set(gca,'xtick',[0:5:60]);
set(gca,'tickdir','out');
colormap(ax4,flipud(cm_Uz));
% colormap(ax4,hot)
% cmocean('-deep');
hold on;
[c,hh]=contour(xc(x_range)*1e-3,zc,1e3*(squeeze(U_z(x_range,1,:,i).^2))',...
    [-0.2:0.1:0.2],'color','k');
set(gca,'ytick',-300:100:0);
set(gca,'xtick',xticks_new);
set(gca,'xticklabel',{[],[],[],[],[],[],[],[],[],'-1','0','1','2','3','4','5'});
% set(gca,'xticklabel',{[round(xc(x_range(1))/1000):round(xc(x_range(end))/1000)]});
z4 = colorbar('location','eastoutside');
set(z4,'position',[0.92 0.046 .02 0.22]); 
hold on; 
fill(xc*1e-3,zb3D,[190 190 190]/225); % topo black
ylim([-300 0]);
xlim([xc(x_range(1)) xc(x_range(end))]*1e-3);
ylabel('Depth (m)','fontsize',18,'fontname','arial')
xxx = xlabel('X (km)','fontsize',18,'fontname','arial');
myticks('x','y',xticks_new(1):.2:xticks_new(end),.8,'popupid');
myticks('y','y',-300:20:0,.8,'popupid');
% set(xxx,'position',get(xxx,'position')+[4 22 1.0]);
set(gca,'fontsize',15);
set(gcf,'color','w');
text(41,-220,'S^2 (10^-^3s^-^2)','fontsize',18,'fontname','arial');

ax5 = axes('position',[0.265 0.769 .15 .08]);
plot((itv(time)-itv(1))/3600,U_lowpass,'k-','linewidth',2);
hold on; plot((itv(time(i))-itv(1))/3600,U_lowpass(i),'r.','markersize',16);
ylim([ 1.0 2.5]);
xlim([0 12])
% set(gca,'xtick',[20:5:35],'xticklabel',{'0','5','10',[]});
ylabel('U (m/s)');xlabel('Time (hour)');grid on;
set(gca,'XAxisLocation','top');
set(gca,'fontsize',10,'xminortick','on','yminortick','on','ticklength',[0.025 0.025]);


frame = getframe(gcf);

writeVideo(temp,frame);
 clf;

end


close(temp);
