clear
load('input/topo_sill.mat','zb3D');
%cd run
ITER = [0:960:225120];
% the model time is 120:120:238200, but deltaT=0.8, so the model iteration should be: 150:150:297750
%ITER = [259200:480:360000]; % 10 minutes of interval: 2 min x 5   [the model only runs for 1366 time steps and it crashes, so 1366*120=163920]
%ITER = [0:480:60000];
%ITER = inf;
%ITER = inf;%itv=inf;
addpath /glade/scratch/liuchihl/temp/MITgcm/utils/matlab
data=rdmnc('run/state.*','X','Y','Z','T','Temp','S','V','U','W',ITER);
%data_phi=rdmnc('run/phiHyd.*','X','Y','Z','T','phiHyd',ITER);

%data=rdmnc('run/3D_diags.*','X','Y','Z','T','THETA','SALT','UVEL','VVEL','WVEL');

%VV=(data.V(:,1:end-1,:,:)+data.V(:,2:end,:,:))/2;
%UU=(data.U(1:end-1,:,:,:)+data.U(2:end,:,:,:))/2;

%[KLdiffKr,ITS]=rdmds('run/KLdiffKr',ITER);
%[KLviscAr,ITS] = rdmds('KLviscAr',ITER);



%imagesc( S.XC, S.YC, S.T(:,:,1)' );

%[U,ITS]=rdmnc('U',itv);
%[V,ITS]=rdmds('V',itv);
%[S,ITS]=rdmds('S',itv);
%[T,ITS]=rdmds('T',itv);
%[W,ITS]=rdmds('W',itv);
%[KLdiffKr,ITS]=rdmds('KLdiffKr',itv);
%[KLviscAr,ITS] = rdmds('KLviscAr',itv);

%cd ..

%save('ustkl_test','T','xc','zb','zc','itv','-v7.3');
%save('Kv','KLdiffKr','-v7.3');

%save('uvt_all','data','zb3D','-v7.3');



xc = data.X; yc = data.Y; zc = data.Z; itv = data.T;
U = data.U(:,1,:,:); V = data.V(:,1,:,:);
W = data.W(:,1,:,:); T = data.Temp(:,1,:,:);
S = data.S(:,1,:,:); %phiHyd = data_phi.phiHyd(:,1,:,:);
save('uvt2D_baro_sponge_wide_orlanski_spinup_ini','U','V','W','T','S','xc','yc','zc','zb3D','itv','-v7.3');
%save('uvt2D_baro_sponge_wide_orlanski_spinup','U','V','W','T','S','xc','yc','zc','zb3D','itv','-v7.3');

%save('uvt2D_baro_sponge_ORlanski_wide_phi','phiHyd','xc','yc','zc','zb3D','itv','-v7.3');

