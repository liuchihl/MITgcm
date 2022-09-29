clear
load('input/topo.mat','xc','zc','zb');
cd run
addpath /home/chihlun/MITgcm/utils/matlab

itv = [172800:120:172920];

%[U,ITS]=rdmds('U',itv);
%[V,ITS]=rdmds('V',itv);
%[S,ITS]=rdmds('S',itv);
%[T,ITS]=rdmds('T',itv);
%[W,ITS]=rdmds('W',itv);
%[KLviscAr,ITS]=rdmds('KLviscAr',itv);
[KLeps,ITS] = rdmds('KLeps',itv);
%[KLdiffKr,ITS]=rdmds('KLdiffKr',itv);

cd ..

%save('ustkl_test','T','xc','zb','zc','itv','-v7.3');
%save('Kv','KLdiffKr','-v7.3');
save('epsilon','KLeps','xc','zb','zc','itv','-v7.3');

