clear all
cd ../teststat2d/debugtestalldom8
tic
load('matlab.mat','XP','xp','T','N','Nx','L','x1','x2','dt')
toc
cd ../..
tic
DD=diffcomp2(XP,T,N,Nx,L,x1,x2,dt);
% DD=diffcomp2(XP,T,N,Nx,L,x1,x2,dt);
toc
cd teststat2d/debugtestalldom8
tic
save('DD.mat','DD','-v7.3')
toc