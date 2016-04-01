clear all
D=0.025;
% x1=4;
x1=4.5;
% x2=6;
x2=5.5;

xm=5;
L=10;
w=1;
T=100;
N=1000000;
Nx=100;
[xp,yp,dt]=stat_2d_diff_romr(D,x1,x2,xm,w,N,T,L,Nx);

x=x1:(x2-x1)/Nx:x2;%close to center
% x=0:L/Nx:L;%all somain
% T=size(xp,2);
t=1:T;
% N=size(xp,1);
C=zeros(Nx,T);
ymin=ones(T,Nx)*w/2;
ymax=ones(T,Nx)*w/2;
for it=1:T
    for i=1:N
        for ix=1:Nx
            if xp(i,it)>=x(ix) && xp(i,it)<x(ix+1)
                C(ix,it)=C(ix,it)+1;
                if ymin(it,ix)>yp(i,it) %&& it==T
                    ymin(it,ix)=yp(i,it);
                end
                if ymax(it,ix)<yp(i,it) %&& it==T
                    ymax(it,ix)=yp(i,it);
                end
            end
        end
    end
end

xx=zeros(1,Nx);
V=zeros(T,Nx);
Vav=zeros(T,Nx);
Cnorm=zeros(Nx,T);

for ix=1:Nx
    xx(ix)=(x(ix)+x(ix+1))/2;
    V(:,ix)=(x(ix+1)-x(ix))*(ymax(:,ix)-ymin(:,ix));
end
for ix=1:Nx
    Vav(:,ix)=(V(:,ix)+V(:,end-ix+1))/2;
    for it=1:T
        Cnorm(ix,it)=C(ix,it)/Vav(it,ix);
    end
end
xxx=zeros(1,Nx-1);
Cnd=zeros(Nx-1,T);
Cd=zeros(Nx-1,T);
for ix=1:Nx-1
    xxx(ix)=(xx(ix)+xx(ix+1))/2;
    for it=1:T
        Cnd(ix,it)=(Cnorm(ix+1,it)-Cnorm(ix,it))/(xx(ix+1)-xx(ix));
        Cd(ix,it)=(C(ix+1,it)-C(ix,it))/(xx(ix+1)-xx(ix));
    end
end



figure(1)
surf(xx,t,Cnorm','Edgecolor','None')
figure(2)
surf(xx,t,C','Edgecolor','None')

% for it=1:T
%     pause(0.1)
%     figure(3)
%     plot(xx,V(it,:))
%     axis([x1 x2 0 (x2-x1)/Nx*w]);
% end
% for it=1:T
%     pause(0.1)
%     figure(4)
%     plot(xx,Vav(it,:))
%     axis([x1 x2 0 (x2-x1)/Nx*w]);
% end
% for it=1:T
%     pause(0.1)
%     figure(5)
%     plot(xx,Cnorm(:,it));
%     axis([x1 x2 0 0.05*T*N]);
% end
% figure(6)
% plot(xx,Cnorm(:,1));
% axis([x1 x2 0 N]);
% figure(7)
% plot(xx,Cnorm(:,end));
% axis([x1 x2 0 N]);

figure(8)
surf(xxx,t,Cnd','Edgecolor','None')
figure(9)
surf(xxx,t,Cd','Edgecolor','None')
figure(12)
kk=0;
for k=0:20:T
    kk=kk+1;
    subplot(3,2,kk)
    if k~=T
        plot(xxx,Cnd(:,k+1)')
    else
        plot(xxx,Cnd(:,T)')
    end
end
figure(11)
kk=0;
for k=0:20:T
    kk=kk+1;
    subplot(3,2,kk)
    if k~=T
        plot(xx,Cnorm(:,k+1)')
    else
        plot(xx,Cnorm(:,T)')
    end
end
