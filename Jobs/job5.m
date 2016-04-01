% tic
cd ..
K=5000;
N=1e4;
Nx=100;
% script_multi_sim_change_dom_trunc
% tic
script_multi_sim_steady_dom_trunc
% toc
cd teststat2d
% mkdir test_multi_change_trunc_less
% cd test_multi_change_trunc_less
% mkdir test_multi_steady_trunc_less2
cd test_multi_steady_trunc_less2
save('matlab.mat','-v7.3')
% toc