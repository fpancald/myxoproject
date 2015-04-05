function [u5,u6]=myxoromr1()
global L T rpol abcdl abcdr DD DE DR r R s1 s1p s2 s3 s4 s4p N0 Nc Div divT icD icE icR abcdlc abcdrc
format long
L = 10;                    %space domain size
xsteps=floor(25*L+1);       %#space points
T=10000;                    %time domain size
tsteps=floor(2*T+1);        %#time points
rpol=1;                     %radius of poles
abcdl=deg3interpol(rpol/2,rpol,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%near left pole
abcdr=deg3interpol(L-rpol,L-rpol/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%near right pole
abcdlc=deg3interpol((L-rpol)/2,L/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%center left 
abcdrc=deg3interpol(L/2,(L+rpol)/2,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%center right
m = 0;                      %cartesian coordinates
x = linspace(0,L,xsteps);   %space domain
t = linspace(0,T,tsteps);   %time domain
%MinCDE-like system parameters
DD=0.28*25/9;
DE=0.6*25/9;
s1=20/9;
s1p=0.028;
s2=0.0063/9;
s3=0.04/9;
s4=0.8/9;
s4p=0.027;
icD=1500;
icE=85;

% RomR parameters
DR=1;%diffusion rate
r=0.1;%association rate
R=0.05;%dissociation rate
N0=400; %basic receptor density 
Nc=1;%molecules per receptor
Div=1;%division ON/OFF switch
divT=T/2;%division start time (end at T)
icR=3000;
%%may need to add change in diffusion at division site

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

u1 = sol(:,:,1);%MinD-like free
u2 = sol(:,:,2);%bounded
u3 = sol(:,:,3);%MinE-like free
u4 = sol(:,:,4);%bounded
u5 = sol(:,:,5);%RomR free
u6 = sol(:,:,6);%bounded

figure(1)
grid off
surf(x,t,u5+u6,'EdgeColor', 'none')
title('RomR')
xlabel('Distance x')
ylabel('Time t')

figure(2)
grid off
surf(x,t,u4,'EdgeColor', 'none')
title('MinE(x)')
xlabel('Distance x')
ylabel('Time t')
% 
% figure(2)
% grid off
% surf(x,t,u3+u4,'EdgeColor', 'none')
% title('MinE(x)')
% xlabel('Distance x')
% ylabel('Time t')
% 
% 
% figure(3);
% for j = 1:100:length(t)
% plot(x,u1(j,:)+u2(j,:),'-k')
% title(['Solution at t =', num2str(t(j))]);
% xlabel('Distance x')
% ylabel('MinD(x)')
% % ylim([-1.1 1.1]);
% 
% pause(0.05);
% end
% 
% figure(4);
% for j = 1:100:length(t)
% plot(x,u3(j,:)+u4(j,:),'-k')
% title(['Solution at t =', num2str(t(j))]);
% xlabel('Distance x')
% ylabel('MinE(x)')
% % ylim([-1.1 1.1]);
% 
% pause(0.05);
% end
% 
% figure(5)
% MinD=u1(1,:)+u2(1,:);
% for j = 2:length(t)
% MinD=MinD+u1(j,:)+u2(j,:);
% end
% % MinD=(u1(1,:)+u2(1,:))/max(u1(1,:)+u2(1,:));
% % for j = 2:length(t)
% % MinD=MinD+(u1(j,:)+u2(j,:))/max(u1(j,:)+u2(j,:));
% % end
% 
% MinD=MinD/length(t);
% plot(x,MinD/max(MinD));
% % plot(x,MinD);
% title('average MinD/maxMinD');
% xlabel('Distance x')
% ylabel('MinD(x)/maxMinD')
% 
% figure(6)
% MinE=u3(1,:)+u4(1,:);
% for j = 2:length(t)
% MinE=MinE+u3(j,:)+u4(j,:);
% end
% % MinE=(u3(1,:)+u4(1,:))/max(u3(1,:)+u4(1,:));
% % for j = 2:length(t)
% % MinE=MinE+(u3(j,:)+u4(j,:))/max(u3(j,:)+u4(j,:));
% % end
% 
% MinE=MinE/length(t);
% plot(x,MinE/max(MinE));
% % plot(x,MinE);
% title('average MinE/maxMinE');
% xlabel('Distance x')
% ylabel('MinE(x)/maxMinE')
% 
% [~,I]=min(MinD/max(MinD));
% % [C,I]=min(MinD);
% x(I)



end

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global L T rpol abcdl abcdr DD DE DR r R s1 s1p s2 s3 s4 s4p N0 Nc Div divT icD icE abcdlc abcdrc

if x<rpol/2 ||x>L-rpol/2 %poles
    
%     N=N0*(1+exp(-u(2)));%exp interaction
    N=N0*(1+u(4)/(4*icE));%linear interaction
    
elseif x>=rpol/2 && x<rpol %near left pole
    
%     ab=deg1interpol(rpol/2,rpol,N0*(1+u(4)/(4*icE)),0);%deg1 receptor spline
%     N=ab(1)*x+ab(2);
    abcd=abcdl*N0*(1+u(4)/(4*icE));
    N=abcd(1)*x^3+abcd(2)*x^2+abcd(3)*x+abcd(4);%deg 3 receptor spline
    
elseif x>L-rpol && x<=L-rpol/2%near right pole
    
%     ab=deg1interpol(L-rpol,L-rpol/2,0,N0*(1+u(4)/(4*icE)));
%     N=ab(1)*x+ab(2);
    abcd=abcdr*N0*(1+u(4)/(4*icE));
    N=abcd(1)*x^3+abcd(2)*x^2+abcd(3)*x+abcd(4);
    
elseif x>(L-rpol)/2 && x<(L+rpol)/2 && Div==1 && t>divT%division ON,dvision started,and near center
    
    if x<L/2%left of center
        
    %     ab=deg1interpol((L-rpol)/2,L/2,0,N0/2*(1+u(4)/(4*icE))*(t-divT)/(T-divT));
    %     N=ab(1)*x+ab(2);
        abcd=abcdlc*N0*(1+u(4)/(4*icE))*(t-divT)/(T-divT);
        N=abcd(1)*x^3+abcd(2)*x^2+abcd(3)*x+abcd(4);
        
    else%right of center
        
    %     ab=deg1interpol((L-rpol)/2,L/2,0,N0/2*(1+u(4)/(4*icE))*(t-divT)/(T-divT));
    %     N=ab(1)*x+ab(2);
        abcd=abcdrc*N0*(1+u(4)/(4*icE))*(t-divT)/(T-divT);
        N=abcd(1)*x^3+abcd(2)*x^2+abcd(3)*x+abcd(4);
    end
    
else%not close to center,not poles or close, division OFF or division not started
    N=0;
end

un=max(N-u(6)/Nc,0)/N0;

c = [1;1;1;1;1;1]; 
f = [DD;1e-10;DE;1e-10;DR;1e-10] .* DuDx; 
s = [-s1*u(1)/(1+s1p*u(4))+s2*u(4)*u(2);+s1*u(1)/(1+s1p*u(4))-s2*u(4)*u(2);-s3*u(1)*u(3)+s4*u(4)/(1+s4p*u(1));+s3*u(1)*u(3)-s4*u(4)/(1+s4p*u(1));-r*un*u(5)+R*u(6);r*un*u(5)-R*u(6)];

end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
global icD icE icR
 %u0=[2*rand*icD;0;2*rand*icE;0;2*rand*icR;0];
 u0=[random('Normal',icD,100);0;random('Normal',icE,25);0;random('Normal',icR,100);0];
 
% --------------------------------------------------------------
end
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = [0;0;0;0;0;0];
ql = [1;1;1;1;1;1];
pr = [0;0;0;0;0;0];
qr = [1;1;1;1;1;1];
end