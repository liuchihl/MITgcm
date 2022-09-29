close all;%figure; 
T(T==0)=nan;
%U(U==0)=nan;
x_range = 750:1377;

%[c,h] = contourf(xc(x_range),zc,squeeze(T(x_range,1,:))',50); set(h,'edgecolor','none');
%colorbar; 
for i=20:40
figure; 
[c,h] = contourf(xc(x_range),zc,squeeze(T(x_range,1,:,i))',50);
%set(h,'edgecolor','none');
colorbar;
%caxis([-1 2]);
%hold on; contourf(xc(x_range),zc,squeeze(U(x_range,1,:,i))',
end

