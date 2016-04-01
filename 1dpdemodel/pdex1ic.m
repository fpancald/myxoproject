function u0 = pdex1ic(x)
global icD icE icR
 %u0=[2*rand*icD;0;2*rand*icE;0;2*rand*icR;0];
 u0=[random('Normal',icD,100);0;random('Normal',icE,25);0;random('Normal',icR,100);0];
 
% --------------------------------------------------------------
end