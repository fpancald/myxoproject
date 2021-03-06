function myxoromrchng(Ntask)
global L T  abcdl  DD DE DR r R s1 s1p s2 s3 s4 s4p N0 Nc Div Diff divT icD icE icR 
rng('shuffle');
seed=rng;
format long
L = 10;                    %max space domain size
xsteps=floor(25*L+1);       %#space points
T=3600*4;%10000                    %time domain size
tsteps=floor(2*T+1);        %#time points

% rpol=1;                     %radius of poles
abcdl=deg3interpol(1/2,1,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%near left pole
% abcdr=deg3interpol(1-rpol,1-rpol/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%near right pole
% abcdlc=deg3interpol((1-rpol)/2,1/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%center left 
% abcdrc=deg3interpol(1/2,(1+rpol)/2,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%center right
m = 0;                      %cartesian coordinates
x = linspace(0,1,xsteps);   %space domain adimensional
t = linspace(0,T,tsteps);   %time domain dimensional
%MinCDE-like system parameters
Ls=6.25;%25
Ts=4.2;%10.8;%8.4;
DD=0.28*Ls/Ts;
DE=0.6*Ls/Ts;
s1=20/Ts;
s1p=0.028;
s2=0.0063/Ts;
s3=0.04/Ts;
s4=0.8/Ts;
s4p=0.027;
icD=1500;
icE=85;

% RomR parameters
DR=1;%diffusion rate
r=0.1;%association rate
R=0.05;%dissociation rate
N0=400; %basic receptor density 
Nc=1;%molecules per receptor
Div=0;%division ON/OFF switch
Diff=0;%spatial diffusion ON/OFF switch
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

filename=strcat(strcat('data',num2str(Ntask)),'.mat');
save(filename)
% figure(1)
% grid off
% surf(x,t,u5+u6,'EdgeColor', 'none')
% title('RomR')
% xlabel('Distance x')
% ylabel('Time t')
% 
% figure(2)
% grid off
% surf(x,t,u4,'EdgeColor', 'none')
% title('MinE(x)')
% xlabel('Distance x')
% ylabel('Time t')

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