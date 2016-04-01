function []=diffusionsave2(Ntask)
tic
rng('shuffle');
scriptarrayjob2
filename=strcat(strcat('DD',num2str(Ntask)),'.mat');
cd teststat2d
mkdir diffusiontestrec
cd diffusiontestrec
save(filename,'DD','-v7.3');
toc
if Ntask==1
    save('data.mat','-v7.3');
end
end