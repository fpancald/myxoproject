function [xp,yp,dt] = stat_2d_diff_romr2(D,x1,x2,xm,w,N,T,L,Nx)
%derive Fick's law-like relation from statistical simulation of proteins
%moving in a 2d domain changing in time and space
% D=0.025;
xl=x1:0.01:xm;
xr=xm:0.01:x2;
y=0:0.01:w;
% w=1;
% x1=4;
% x2=6;
% xm=5;
dx=(x2-x1)/Nx;
dt=dx^2/(2*D);
sigma=dx;
%poles shape
B=w/2;
R=w/2;
Al=w/2;
xlp=0:0.01:Al;
Slpu(1,:)=(sqrt(R^2-(xlp-Al).^2)+B);
Slpd(1,:)=1-(sqrt(R^2-(xlp-Al).^2)+B);
Ar=L-w/2;
xrp=Ar:0.01:L;
Srpu(1,:)=(sqrt(R^2-(xrp-Ar).^2)+B);
Srpd(1,:)=1-(sqrt(R^2-(xrp-Ar).^2)+B);

% ym=ones(1,T+1)*w/2;

% %parabolic profile
% ym=w:-(w/2)/T:w/2;
% for i=1:T
%     b(i)=(1+x1-xm-ym(i+1)^2)/(2*(1-ym(i+1)));
%     a(i)=(ym(i+1)-b(i))^2+xm;
%     
%     Slu(i,:)=sqrt(a(i)-xl)+b(i);
%     Sld(i,:)=1-(sqrt(a(i)-xl)+b(i));
%     Sru(i,:)=sqrt(a(i)-(2*xm-xr))+b(i);
%     Srd(i,:)=1-(sqrt(a(i)-(2*xm-xr))+b(i));
% end

%circular profile
r=100*w/2*ones(1,T+1);
% r=10*w/2:-(9*w/2)/T:w/2;
for t=1:T
    b(t)=w/2;
    a(t)=x1;
%     r(t+1)
    ynorm(t)=sqrt(r(t+1)^2-(x1-a(t)).^2)+b(t);
    Slu(t,:)=(sqrt(r(t+1)^2-(xl-a(t)).^2)+b(t))/ynorm(t);
    Sld(t,:)=1-(sqrt(r(t+1)^2-(xl-a(t)).^2)+b(t))/ynorm(t);
    Sru(t,:)=(sqrt(r(t+1)^2-((2*xm-xr)-a(t)).^2)+b(t))/ynorm(t);
    Srd(t,:)=1-(sqrt(r(t+1)^2-((2*xm-xr)-a(t)).^2)+b(t))/ynorm(t);
end

% b=0.5;
% for i=1:length(ym)
%     a=(x1^2-xm^2-ym(i)^2+2*ym(i)*b)/(2*(x1-xm));
%     r=sqrt((x1-a)^2+b^2);
%     
%     Slu(i,:)=sqrt(r^2-(xl-a).^2)+b;
%     Sld(i,:)=-sqrt(r^2-(xl-a).^2)+b;
%     
%     a=(x2^2-xm^2-ym(i)^2+2*ym(i)*b)/(2*(x2-xm));
%     r=sqrt((x2-a)^2+b^2);
%     Sru(i,:)=sqrt(r^2-(xr-a).^2)+b;
%     Srd(i,:)=-sqrt(r^2-(xr-a).^2)+b;
% end

% r=0.5:0.1:1;
% for ir=1:length(r)
%     Slu(ir,:)=sqrt(r(ir)^2-(xl-x1).^2)+w/2;
%     Sld(ir,:)=-sqrt(r(ir)^2-(xl-x1).^2)+w/2;
%     Sru(ir,:)=sqrt(r(ir)^2-(xr-x2).^2)+w/2;
%     Srd(ir,:)=-sqrt(r(ir)^2-(xr-x2).^2)+w/2;
% end

figure(10)
hold on
% plot(xl,Slu(end,:))
% plot(xl,Sld(end,:))
% plot(xr,Sru(end,:))
% plot(xr,Srd(end,:))
for t=1:T
    plot(xl,Slu(t,:))
    plot(xl,Sld(t,:))
    plot(xr,Sru(t,:))
    plot(xr,Srd(t,:))
    if t==1
        plot(xlp,Slpu(t,:))
        plot(xlp,Slpd(t,:))
        plot(xrp,Srpu(t,:))
        plot(xrp,Srpd(t,:))
        xc=Al:0.01:Ar;
        plot(xc,w*ones(size(xc)))
        plot(xc,zeros(size(xc)))
        axis([-w/2 L+w/2 -w/2 w+w/2])%all domain
%         axis([x1 x2 -w/2 w+w/2])%close to center
        daspect([1 1 1])
    end
end
hold off
% xp(1)=2*rand+4;
% yp(1)=rand;

% N=1000;
% T=10000;
xp=zeros(N,T);
yp=zeros(N,T);
for n=1:N
    t=1;
while t<T+1
    if t==1;
        %all domain
        xnew=rand*L;
        ynew=rand*w;

%         %close to center
%         xnew=rand*(x2-x1)+x1;
%         ynew=rand*w;
        
    else
    
%         angle=rand*pi;
%         dist=randn*D;
% 
%         xnew=xp(n,t-1)+cos(angle)*dist;
%         ynew=yp(n,t-1)+sin(angle)*dist;
        
%         distx=randn*D/sqrt(2);
%         disty=randn*D/sqrt(2);
        
        distx=randn*sigma;
        disty=randn*sigma;


        xnew=xp(n,t-1)+distx;
        ynew=yp(n,t-1)+disty;
    end
    
%     if xnew<4
%         xnew=6-(4-xnew);
%     elseif xnew>6
%         xnew=4+(xnew-6);
%     end
%     xp(i,t)=xnew;
%     yp(i,t)=ynew;

%all domain
if xnew<0
        t=t-1;
    elseif xnew>L
        t=t-1;

% %close to center
% if xnew<x1
%         t=t-1;
%     elseif xnew>x2
%         t=t-1;
else
    xp(n,t)=xnew;
    yp(n,t)=ynew;
    
    if ynew<0 || ynew>w
        t=t-1;
    elseif ynew>w/2 && xnew>x1 && xnew<xm
%         ylim=sqrt(a(t)-xnew)+b(t);%parabola
        ylim=(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);%circle
        if ynew>ylim
            t=t-1;
        end
    elseif ynew<w/2 && xnew>x1 && xnew<xm
%         ylim=1-(sqrt(a(t)-xnew)+b(t));%parabola
        ylim=1-(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);%circle
        if ynew<ylim
            t=t-1;
        end
    elseif ynew>w/2 && xnew>xm && xnew<x2
%         ylim=sqrt(a(t)-(2*xm-xnew))+b(t);%parabola
        ylim=(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t);%circle
        if ynew>ylim
            t=t-1;
        end
    elseif ynew<w/2 && xnew>xm && xnew<x2
%         ylim=1-(sqrt(a(t)-(2*xm-xnew))+b(t));%parabola
        ylim=1-(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t);%circle
        if ynew<ylim
            t=t-1;
        end
        %uncomment if all domain
    elseif ynew>w/2 && xnew<Al
        ylim=(sqrt(R^2-(xnew-Al).^2)+B);%circle
        if ynew<ylim
            t=t-1;
        end
    elseif ynew<w/2 && xnew<Al
        ylim=1-(sqrt(R^2-(xnew-Al).^2)+B);%circle
        if ynew<ylim
            t=t-1;
        end
    elseif ynew>w/2 && xnew<Ar
        ylim=(sqrt(R^2-(xnew-Ar).^2)+B);%circle
        if ynew<ylim
            t=t-1;
        end
    elseif ynew<w/2 && xnew<Ar
        ylim=1-(sqrt(R^2-(xnew-Ar).^2)+B);%circle
        if ynew<ylim
            t=t-1;
        end
    end
    
end
    
    t=t+1;
end
end
% figure()
% plot(xp,yp)
% hold on
% for i=1:ceil(N/10):N
%     figure(i)
%     hold on

%     plot(xl,Slu(end,:))
%     plot(xl,Sld(end,:))
%     plot(xr,Sru(end,:))
%     plot(xr,Srd(end,:))

%     plot(xp(i,:),yp(i,:))
%     scatter(xp(i,:),yp(i,:))
    
%     scatter(xp(i,1),yp(i,1));

%     hold off

% end

% hold off