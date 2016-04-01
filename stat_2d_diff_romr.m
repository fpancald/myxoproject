%Derive Fick's law-like relation from statistical simulation of proteins
%moving in a 2d domain changing in time and space
%use repulsive condition at the boundary (i.e. if the step goes out of the
%boundary re/select the step)
%Input:
%D:diffusion coefficient in 2D domain
%x1,x2: left and right bound for the zone close to the center (x-coord.)
%xm: center (x-coord.)
%N: number of particles to simulate
%T: number of time steps to simulate
%L: total domain length (x-coord.)
%Nx: number of intervals in the x direction (around the center)
%state:0 unchanged, 1 (or not zero) changing
%Output:
%xp,yp: NxT matrices containing the x and y coord. of each particle at each time step
%dt: time step length correspondent to the diffusion coeffcient and lenght of space step dx chosen
function [xp,yp,dt] = stat_2d_diff_romr(D,x1,x2,xm,w,N,T,L,Nx,state)
xl=x1:0.01:xm;  %left of center x-coord.
xr=xm:0.01:x2;  %right of center x-coord
y=0:0.01:w;     %y domain 
dx=(x2-x1)/Nx;  %space step
dt=dx^2/(2*D);  %time step
sigma=dx;       %std of random step in 1 dimension

%poles shape (set the equations for poles at domain extremities)
%not yet used in simulation
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

% %parabolic profile
% ym=ones(1,T+1)*w/2;

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

%circular profile (ellipsis in x directions become semicircles at final time)
if state==0
    r=100*w/2*ones(1,T+1);% constant and (almost rectangular)unchanged domain
else
    r=10*w/2:-(9*w/2)/T:w/2; %changing domain
end
b=zeros(1,T);
a=zeros(1,T);
ynorm=zeros(1,T);
Slu=zeros(T,length(xl));
Sld=zeros(T,length(xl));
Sru=zeros(T,length(xr));
Srd=zeros(T,length(xr));
for t=1:T
    b(t)=w/2;
    a(t)=x1;
    ynorm(t)=sqrt(r(t+1)^2-(x1-a(t)).^2)+b(t);
    Slu(t,:)=(sqrt(r(t+1)^2-(xl-a(t)).^2)+b(t))/ynorm(t);
    Sld(t,:)=w-(sqrt(r(t+1)^2-(xl-a(t)).^2)+b(t))/ynorm(t);
    Sru(t,:)=(sqrt(r(t+1)^2-((2*xm-xr)-a(t)).^2)+b(t))/ynorm(t);
    Srd(t,:)=w-(sqrt(r(t+1)^2-((2*xm-xr)-a(t)).^2)+b(t))/ynorm(t);
end

%plot domain
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
%         axis([-w/2 L+w/2 -w/2 w+w/2])%all domain
        axis([x1 x2 -w/2 w+w/2])%close to center
        daspect([1 1 1])
    end
end
hold off

%simulation of random steps
xp=zeros(N,T);
yp=zeros(N,T);
for n=1:N
    t=1;
while t<T+1
    if t==1;
%         %all domain
%         xnew=rand*L;
%         ynew=rand*w;

        %close to center
        xnew=rand*(x2-x1)+x1;
        ynew=rand*w;
        
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
    
%check for new position inside current boundaries
%1st check x-coord left and right bounds
%2nd check y-coord overall upper and lower bound
%3rd check (for each sector that local ylimit is not excedeed (upper or lower)

% %all domain
% if xnew<0
%         t=t-1;
%     elseif xnew>L
%         t=t-1;

%close to center
if xnew<x1
        t=t-1;
    elseif xnew>x2
        t=t-1;
        
        
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
        ylim=w-(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);%circle
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
        ylim=w-(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t);%circle
        if ynew<ylim
            t=t-1;
        end
        
        %uncomment if all domain
%     elseif ynew>w/2 && xnew<Al
%         ylim=(sqrt(R^2-(xnew-Al).^2)+B);%circle
%         if ynew<ylim
%             t=t-1;
%         end
%     elseif ynew<w/2 && xnew<Al
%         ylim=1-(sqrt(R^2-(xnew-Al).^2)+B);%circle
%         if ynew<ylim
%             t=t-1;
%         end
%     elseif ynew>w/2 && xnew<Ar
%         ylim=(sqrt(R^2-(xnew-Ar).^2)+B);%circle
%         if ynew<ylim
%             t=t-1;
%         end
%     elseif ynew<w/2 && xnew<Ar
%         ylim=1-(sqrt(R^2-(xnew-Ar).^2)+B);%circle
%         if ynew<ylim
%             t=t-1;
%         end
    end
    
end
    
    t=t+1;
end
end

%plot particles position
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