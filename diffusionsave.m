function []=diffusionsave(Ntask)
tic
rng('shuffle');
scriptarrayjob
filename=strcat(strcat('DD',num2str(Ntask)),'.mat');
cd teststat2d
mkdir diffusiontest2
cd diffusiontest2
save(filename,'DD','-v7.3');
toc
end