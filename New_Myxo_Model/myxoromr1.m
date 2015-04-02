function [u5,u6]=myxoromr1()
global L rpol abcdl abcdr
format long
L = 10;
rpol=1;
abcdl=deg3interpol(rpol/2,rpol,1,0);
abcdr=deg3interpol(L-rpol,L-rpol/2,0,1);
m = 0;
x = linspace(0,L,500);
t = linspace(0,5000,10000);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);
u5 = sol(:,:,5);
u6 = sol(:,:,6);

figure(1)
grid off
surf(x,t,u5+u6,'EdgeColor', 'none')
title('RomR')
xlabel('Distance x')
ylabel('Time t')

figure(2)
grid off
surf(x,t,u4,'EdgeColor', 'none')
title('MinD(x)')
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
global L rpol abcdl abcdr
DD=0.28*25/9;
DE=0.6*25/9;
DR=1;
r=0.1;
R=0.05;
s1=20/9;
s1p=0.028;
s2=0.0063/9;
s3=0.04/9;
s4=0.8/9;
s4p=0.027;

if x<rpol/2 ||x>L-rpol/2
%     N=400*(1+exp(-u(2)));
    N=400*(1+u(4)/340);
elseif x>=rpol/2 && x<rpol
%     ab=deg1interpol(rpol/2,rpol,400*(1+u(4)/340),0);
%     N=ab(1)*x+ab(2);
    abcd=abcdl*400*(1+u(4)/340);
    N=abcd(1)*x^3+abcd(2)*x^2+abcd(3)*x+abcd(4);
elseif x>L-rpol && x<=L-rpol/2
%     ab=deg1interpol(L-rpol,L-rpol/2,0,400*(1+u(4)/340));
%     N=ab(1)*x+ab(2);
    abcd=abcdr*400*(1+u(4)/340);
    N=abcd(1)*x^3+abcd(2)*x^2+abcd(3)*x+abcd(4);
else 
    N=0;
end
Nc=1;
un=max(N-u(6)/Nc,0)/400;

c = [1;1;1;1;1;1]; 
f = [DD;1e-10;DE;1e-10;DR;1e-10] .* DuDx; 
s = [-s1*u(1)/(1+s1p*u(4))+s2*u(4)*u(2);+s1*u(1)/(1+s1p*u(4))-s2*u(4)*u(2);-s3*u(1)*u(3)+s4*u(4)/(1+s4p*u(1));+s3*u(1)*u(3)-s4*u(4)/(1+s4p*u(1));-r*un*u(5)+R*u(6);r*un*u(5)-R*u(6)];

end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
 %u0=[2*rand*1500;0;2*rand*85;0;2*rand*3000;0];
 u0=[random('Normal',1500,100);0;random('Normal',85,25);0;random('Normal',3000,100);0];
 
% --------------------------------------------------------------
end
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = [0;0;0;0;0;0];
ql = [1;1;1;1;1;1];
pr = [0;0;0;0;0;0];
qr = [1;1;1;1;1;1];
end