clear all
cd ../teststat2d/debugtestalldom7
tic
load('matlab.mat','XP','xp','T','N','Nx','L','x1','x2','dt')
toc
cd ../..
tic
DD=diffcomp2(XP,T,N,Nx*4,L,x1,x2,dt);
% DD=diffcomp2(XP,T,N,Nx,L,x1,x2,dt);
toc
cd teststat2d/debugtestalldom7
tic
save('DD_finer.mat','DD','-v7.3')
toc
