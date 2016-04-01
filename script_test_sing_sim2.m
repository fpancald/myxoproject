
%single simulation test script
clear all
tic
D=0.025;
x1=4.5;
x2=5.5;
xm=5;
L=10;
w=1;
T=100;
N=10000;
Nx=100;
K=10;
% x=x1:(x2-x1)/Nx:x2;%close to center
x=0:(x2-x1)/Nx:L;%all domain

state=1;
for k=1:K
    [xp,yp,dt]=stat_2d_diff_romr6(D,x1,x2,xm,w,N,T,L,Nx,state);
    XP{k}=xp;
end

cd teststat2d
mkdir debugtestalldom9
cd debugtestalldom9
save('matlab.mat','-v7.3')
toc
cd ..
cd ..
