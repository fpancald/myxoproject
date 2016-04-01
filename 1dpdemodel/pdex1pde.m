% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global T L DD DE DR N0  Nc icE s1 s2 s3 s4 s1p s4p r R abcdl
    RL=(L/2*(1+t/T));
    RLsq=RL^2; %current domain length squared
    DDtx=DD/RLsq;   %diff coeff in reference domain
    DEtx=DE/RLsq;
    DRtx=DR/RLsq;
    DDtxDx=0;
    DEtxDx=0;
    DRtxDx=0;
    rpol=1/RL;


if x<rpol/2 ||x>1-rpol/2 %poles
    
%     N=N0*(1+exp(-u(2)));%exp interaction
    N=N0*(1+u(4)/(4*icE));%linear interaction
    
elseif x>=rpol/2 && x<rpol %near left pole
    
%     ab=deg1interpol(rpol/2,rpol,N0*(1+u(4)/(4*icE)),0);%deg1 receptor spline
%     N=ab(1)*x+ab(2);
    abcd=abcdl*N0*(1+u(4)/(4*icE));
    N=abcd(1)*(x*RL)^3+abcd(2)*(x*RL)^2+abcd(3)*(x*RL)+abcd(4);%deg 3 receptor spline
    
elseif x>1-rpol && x<=1-rpol/2%near right pole
    
%     ab=deg1interpol(L-rpol,L-rpol/2,0,N0*(1+u(4)/(4*icE)));
%     N=ab(1)*x+ab(2);
    abcd=abcdl*N0*(1+u(4)/(4*icE));
    N=abcd(1)*((1-x)*RL)^3+abcd(2)*((1-x)*RL)^2+abcd(3)*((1-x)*RL)+abcd(4);
    
else%not close to center,not poles or close, division OFF or division not started
    N=0;
end

un=max(N-u(6)/Nc,0)/N0;

c = [1;1;1;1;1;1]; 
f = [DDtx;1e-10;DEtx;1e-10;DRtx;1e-10] .* DuDx;%+[DDtxDx;0;DEtxDx;0;DRtxDx;0] .*u; 
s = [-s1*u(1)/(1+s1p*u(4))+s2*u(4)*u(2);+s1*u(1)/(1+s1p*u(4))-s2*u(4)*u(2);-s3*u(1)*u(3)+s4*u(4)/(1+s4p*u(1));+s3*u(1)*u(3)-s4*u(4)/(1+s4p*u(1));-r*un*u(5)+R*u(6);r*un*u(5)-R*u(6)]+[DDtxDx;0;DEtxDx;0;DRtxDx;0] .*DuDx;

end
% --------------------------------------------------------------