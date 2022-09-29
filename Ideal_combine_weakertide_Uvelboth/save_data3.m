clear
load('input/topo_sill.mat','xc','yc','zc','zb');
cd run
addpath /glade/scratch/liuchihl/temp/MITgcm/utils/matlab

%addpath ../prog/
%itv=inf;
%itv = 64800;
%itv=[18000:240:26400];
%itv=[86400:120:237600];
%%%itv = [144000:120:302400];
%itv = [228120:40:235320];
%itv = [204960:40:212160];
itv = [216000:400:360000];   % start from 1.5 to 2.5 days
%itv=inf;
[U,ITS]=rdmds('U',itv);
%[V,ITS]=rdmds('V',itv);
%[S,ITS]=rdmds('S',itv);
[T,ITS]=rdmds('T',itv);
[W,ITS]=rdmds('W',itv);
[KLviscAr,ITS]=rdmds('KLviscAr',itv);

%[KLdiffKr,ITS]=rdmds('KLdiffKr',itv);

%U = U(:,1,:,1:1321);
%S = S(:,1,:,1:1321);
%T = T(:,1,:,1:1321);
%W = W(:,1,:,1:1321);
%KLviscAr = KLviscAr(:,1,:,1:1321);
%KLdiffKr = KLdiffKr(:,1,:,1:1321);


%[Eta,ITS]=rdmds('Eta',itv);
cd ..

save('uwtAz_combine_weakertide_bothU','U','T','W','KLviscAr','xc','zb','zc','ITS','-v7.3');
%save('Kv','KLdiffKr','-v7.3');
%save('ustkl','U','W','T','S','KLviscAr','xc','zb3D','zc','itv','-v7.3');
%save('LOW_TIDE','U','W','T','S','KLviscAr','xc','zb','zc','itv_l','-v7.3');

