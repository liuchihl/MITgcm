load('uvt_all.mat');
xc = data.X; yc = data.Y; zc = data.Z; itv = data.T;
U = data.U(:,length(yc)/2,:,:); V = data.V(:,length(yc)/2,:,:);
W = data.W(:,length(yc)/2,:,:); T = data.Temp(:,length(yc)/2,:,:);
clear data
save('uvt2D','U','V','W','T','xc','yc','zc','itv','-v7.3');

