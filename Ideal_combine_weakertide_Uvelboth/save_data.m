clear
%%%% mpowr !
cd run
%addpath('/home/chihlun/MITgcm_c65r/YILAN ridge/advect_xz_nonhydro/input');
addpath /home/chihlun/MITgcm/utils/matlab
addpath /home/chihlun/MITgcm/mhlib/mhlib/mastermatlab
[U,ITS]=rdmds('U',NaN);
[V,ITS]=rdmds('V',NaN);
[W,ITS]=rdmds('W',NaN);
[T,ITS]=rdmds('T',NaN);
[KLviscAr,ITS]=rdmds('KLviscAr',NaN);
%[Eta,ITS]=rdmds('Eta',NaN);
cd ..

load('input/topo.mat','xc','zc');
%X=cumsum(delX);
%Y=cumsum(delY)';
%Z=cumsum(delZ)';


%dim=size(U);

%S1=nan(dim(1),dim(2),dim(3),dim(4));
%S2=nan(dim(1),dim(2),dim(3),dim(4));
%for i=1:dim(4)  %time
%i
%for j=1:dim(1)  % X-axis
%ut=squeeze(U(j,:,:,i));
%vt=squeeze(V(j,:,:,i));
%S1(j,:,:,i)=mmderiv(-Z,ut')';
%S2(j,:,:,i)=mmderiv(-Z,vt')';
%end
%end

%shear=log10(S1.^2+S2.^2);

%vort=nan(dim(1),dim(2),dim(3),dim(4));
%divg=nan(dim(1),dim(2),dim(3),dim(4));
%for i=1:dim(4)  %time
%i
%for j=1:dim(3)  %depth
%ut=squeeze(U(:,:,j,i));
%vt=squeeze(V(:,:,j,i));
%vort(:,:,j,i)=mmderiv(X,vt)-mmderiv(Y',ut')';
%divg(:,:,j,i)=mmderiv(X,ut)+mmderiv(Y',vt')';
%end
%end
save('UVWT','T','U','W','KLviscAr','ITS','xc','zc','-v7.3');
%save('UVWT','T','lon','lat','ITS','X','Y','Z','-v7.3');

