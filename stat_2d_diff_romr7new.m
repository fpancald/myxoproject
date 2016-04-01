%added receptors 1st instance
%Derive Fick's law-like relation from statistical simulation of proteins
%moving in a 2d domain changing in time and space
%use Reflective condition at the boundary (i.e. incidence angle respect to normal to surface is specular to reflective angle)
%Whole  cell circular profile
%Updated I.C.
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
function [xp,yp,dt] = stat_2d_diff_romr7new(D,x1,x2,xm,w,N,T,L,Nx,state,Type,Rswitch)

%uncomment for plot
% xl=x1:0.01:xm;  %left of center x-coord.
% xr=xm:0.01:x2;  %right of center x-coord
% y=0:0.01:w;     %y domain 

dx=(x2-x1)/Nx;  %space step
dt=dx^2/(2*D);  %time step                          %should we bind together space and time step through usual diffusion relation?
sigma=dx;       %std of random step in 1 dimension
xtol=dx/1000;

%poles shape (set the equations for poles at domain extremities)

B=w/2;
R=w/2;
Al=w/2;

%uncomment for plot
% xlp=0:0.01:Al;
% Slpu(1,:)=(sqrt(R^2-(xlp-Al).^2)+B);
% Slpd(1,:)=w-(sqrt(R^2-(xlp-Al).^2)+B);

Ar=L-w/2;

%uncomment for plot
% xrp=Ar:0.01:L;
% Srpu(1,:)=(sqrt(R^2-(xrp-Ar).^2)+B);
% Srpd(1,:)=w-(sqrt(R^2-(xrp-Ar).^2)+B);

%circular profile (ellipsis in x directions become semicircles at final time)
if state==0
    r=100*w/2*ones(1,T+1);% constant and (almost rectangular)unchanged domain
else
    r=10*w/2:-(9*w/2)/T:w/2; %changing domain
end
b=zeros(1,T);
a=zeros(1,T);
ynorm=zeros(1,T);

%uncomment for plot 
% Slu=zeros(T,length(xl));
% Sld=zeros(T,length(xl));
% Sru=zeros(T,length(xr));
% Srd=zeros(T,length(xr));

for t=1:T
    b(t)=w/2;
    a(t)=x1;
    ynorm(t)=sqrt(r(t+1)^2-(x1-a(t)).^2)+b(t);
    
    %uncomment for plot 
%     Slu(t,:)=(sqrt(r(t+1)^2-(xl-a(t)).^2)+b(t))/ynorm(t);
%     Sld(t,:)=w-(sqrt(r(t+1)^2-(xl-a(t)).^2)+b(t))/ynorm(t);
%     Sru(t,:)=(sqrt(r(t+1)^2-((2*xm-xr)-a(t)).^2)+b(t))/ynorm(t);
%     Srd(t,:)=w-(sqrt(r(t+1)^2-((2*xm-xr)-a(t)).^2)+b(t))/ynorm(t);
    
end

%uncomment for plot 
% %plot domain
% figure(10)
% hold on
% % plot(xl,Slu(end,:))
% % plot(xl,Sld(end,:))
% % plot(xr,Sru(end,:))
% % plot(xr,Srd(end,:))
% for t=1:T
%     plot(xl,Slu(t,:))
%     plot(xl,Sld(t,:))
%     plot(xr,Sru(t,:))
%     plot(xr,Srd(t,:))
%     if t==1
%         plot(xlp,Slpu(t,:))
%         plot(xlp,Slpd(t,:))
%         plot(xrp,Srpu(t,:))
%         plot(xrp,Srpd(t,:))
%         xc=Al:0.01:Ar;
%         plot(xc,w*ones(size(xc)))
%         plot(xc,zeros(size(xc)))
% %         axis([-w/2 L+w/2 -w/2 w+w/2])%all domain
%         axis([x1 x2 -w/2 w+w/2])%close to center
%         daspect([1 1 1])
%     end
% end
% hold off

%receptors distribution
if Rswitch~=0
	perc=0.165;
	[Lr,Rr,rstep]=receptors(perc,N,w);
	%bounded state
	bound=zeros(N); %0 particle not bound 1 particle bound
	%binding/unbinding probabilities
	bp=0.2;
	up=0.0001;
end
%simulation of random steps
xp=zeros(N,T);
yp=zeros(N,T);
for n=1:N
    t=1;
while t<T+1
    if t==1 
        %I.C.
	if Type==0
		 xnew=rand*L;
		 ynew=rand*w;
	elseif Type==1
		xnew=rand*L/2;
		ynew=rand*w;
	elseif Type==2
		if rand<0.6
			xnew=rand*L/2;
			ynew=rand*w;
		else
			xnew=rand*L/2+L/2;
			ynew=rand*w;
		end
	end

         if xnew<Al
             ymax=(sqrt(R^2-(xnew-Al).^2)+B);
             ymin=w-(sqrt(R^2-(xnew-Al).^2)+B);
         elseif xnew>Ar
             ymax=(sqrt(R^2-(xnew-Ar).^2)+B);
             ymin=w-(sqrt(R^2-(xnew-Ar).^2)+B);
         else
             ymin=0;
             ymax=w;
         end
         if ynew<=ymin || ynew>=ymax
             t=0;
         end
         
         %imaginary check
%          if imag(xnew)~=0 || imag(ynew)~=0
%              xnew
%              ynew
%              t
%              pause
%          end


        %end I.C.
        
    else
    
     %check that old position is still inside new domain boundary. If not
     %the y coordinate is corrected accordingly (different approaches are
     %possible here)
        xold=xp(n,t-1);
        yold=yp(n,t-1);
        
        if xold<=0 %need to think of better fix %this can happen because of the tollerance in refl and trunc
            xold=xtol;
        elseif xold>=L
            xold=L-xtol;
        end
        
        if xold<Al
             ymax=(sqrt(R^2-(xold-Al).^2)+B);
             ymin=w-(sqrt(R^2-(xold-Al).^2)+B);
         elseif xold>Ar
             ymax=(sqrt(R^2-(xold-Ar).^2)+B);
             ymin=w-(sqrt(R^2-(xold-Ar).^2)+B);
        elseif xold>x1 && xold<=xm
            ymax=(sqrt(r(t+1)^2-(xold-a(t)).^2)+b(t))/ynorm(t);
            ymin=w-(sqrt(r(t+1)^2-(xold-a(t)).^2)+b(t))/ynorm(t);
        elseif xold<x2 && xold>xm
            ymax=(sqrt(r(t+1)^2-((2*xm-xold)-a(t)).^2)+b(t))/ynorm(t);
            ymin=w-(sqrt(r(t+1)^2-((2*xm-xold)-a(t)).^2)+b(t))/ynorm(t); 
        else
             ymin=0;
             ymax=w;
        end
        
        if yold>=ymax
            yold=ymax-xtol;
        elseif yold<=ymin
            yold=ymin+xtol;
        end
        
        if Rswitch~=0
            %binding and release
            if xold<w/2
                [ip,jp]=locator2d(xold,0,w/2,rstep,yold,0,w,rstep);
                if Lr(ip,jp)>0 && bound(n)==0
                    if rand<bp
                        bound(n)=1;
                        Lr(ip,jp)=Lr(ip,jp)-1;
                    end
                elseif bound(n)==1;
                    if rand<up
                        bound(n)=0;
                        Lr(ip,jp)=Lr(ip,jp)+1;
                    end
                end
            elseif xold>L-w/2
                [ip,jp]=locator2d(xold,L-w/2,L,rstep,yold,0,w,rstep);
                if Rr(ip,jp)>0 && bound(n)==0
                    if rand<bp
                        bound(n)=1;
                        Rr(ip,jp)=Rr(ip,jp)-1;
                    end
                elseif bound(n)==1;
                    if rand<up
                        bound(n)=0;
                        Rr(ip,jp)=Rr(ip,jp)+1;
                    end
                end
            end


            if bound(n)==0
            %         distx=randn*D/sqrt(2);
            %         disty=randn*D/sqrt(2);
                distx=randn*sigma;
                disty=randn*sigma;
                xnew=xold+distx;
                ynew=yold+disty;
            else
                xnew=xold;
                ynew=yold;
            end
        else
            distx=randn*sigma;
            disty=randn*sigma;
            xnew=xold+distx;
            ynew=yold+disty;
        end


        %immaginary check
    %         if imag(xnew)~=0 || imag(ynew)~=0
    %             0
    %              xnew
    %              ynew
    %              t
    %              xold
    %              yold
    %              xp(n,t-2)
    %              yp(n,t-2)
    %              distx
    %              disty
    %              xp(n,1)
    %              yp(n,1)
    %              pause
    %          end

            damp=1;%linear damping coefficient for reflection length
    %         [xnew,ynew]=reflectbc(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %check for new position inside current boundaries
    %1st check x-coord left and right bounds
    %2nd check y-coord overall upper and lower bound
    %3rd check (for each sector that local ylimit is not excedeed (upper or lower)

            if xnew<=0 %need to think of better fix %this can happen because of the tollerance in refl and trunc
                [xnew,ynew]=reflectbc3(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
            elseif xnew>=L
                [xnew,ynew]=reflectbc3(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
            elseif (xnew>xm && xold<xm) || (xnew<xm && xold>xm) %this is to make sure there is no intersection with boundary passing from one side to the other
                [xnew,ynew]=reflectbc3(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
            end

            if xnew<Al
                 ymax=(sqrt(R^2-(xnew-Al).^2)+B);
                 ymin=w-(sqrt(R^2-(xnew-Al).^2)+B);
             elseif xnew>Ar
                 ymax=(sqrt(R^2-(xnew-Ar).^2)+B);
                 ymin=w-(sqrt(R^2-(xnew-Ar).^2)+B);
            elseif xnew>x1 && xnew<=xm
                ymax=(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);
                ymin=w-(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);
            elseif xnew<x2 && xnew>xm
                ymax=(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t);
                ymin=w-(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t); 
            else
                 ymin=0;
                 ymax=w;
            end

            if ynew>=ymax
                [xnew,ynew]=reflectbc3(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
            elseif yold<=ymin
                [xnew,ynew]=reflectbc3(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
            end

    %         [xnew,ynew]=reflectbc3(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %         tol=1e-6;
    %         if xnew<x1-tol || xnew>x2+tol
    %                 [xnew,ynew]=reflectbc2(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %         else
    %             if ynew<0 || ynew>w
    %                 [xnew,ynew]=reflectbc2(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %             elseif ynew>w/2 && xnew>x1 && xnew<xm
    %         %         ylim=sqrt(a(t)-xnew)+b(t);%parabola
    %                 ylim=(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);%circle
    %                 if ynew>ylim
    %                     [xnew,ynew]=reflectbc2(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %                 end
    %             elseif ynew<w/2 && xnew>x1 && xnew<xm
    %         %         ylim=1-(sqrt(a(t)-xnew)+b(t));%parabola
    %                 ylim=w-(sqrt(r(t+1)^2-(xnew-a(t)).^2)+b(t))/ynorm(t);%circle
    %                 if ynew<ylim
    %                     [xnew,ynew]=reflectbc2(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %                 end
    %             elseif ynew>w/2 && xnew>xm && xnew<x2
    %         %         ylim=sqrt(a(t)-(2*xm-xnew))+b(t);%parabola
    %                 ylim=(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t);%circle
    %                 if ynew>ylim
    %                     [xnew,ynew]=reflectbc2(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %                 end
    %             elseif ynew<w/2 && xnew>xm && xnew<x2
    %         %         ylim=1-(sqrt(a(t)-(2*xm-xnew))+b(t));%parabola
    %                 ylim=w-(sqrt(r(t+1)^2-((2*xm-xnew)-a(t)).^2)+b(t))/ynorm(t);%circle
    %                 if ynew<ylim
    %                     [xnew,ynew]=reflectbc2(a(t),b(t),r(t+1),ynorm(t),x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp);
    %                 end
    %                 
    %                 %uncomment if all domain
    %         %     elseif ynew>w/2 && xnew<Al
    %         %         ylim=(sqrt(R^2-(xnew-Al).^2)+B);%circle
    %         %         if ynew<ylim
    %         %             t=t-1;
    %         %         end
    %         %     elseif ynew<w/2 && xnew<Al
    %         %         ylim=1-(sqrt(R^2-(xnew-Al).^2)+B);%circle
    %         %         if ynew<ylim
    %         %             t=t-1;
    %         %         end
    %         %     elseif ynew>w/2 && xnew<Ar
    %         %         ylim=(sqrt(R^2-(xnew-Ar).^2)+B);%circle
    %         %         if ynew<ylim
    %         %             t=t-1;
    %         %         end
    %         %     elseif ynew<w/2 && xnew<Ar
    %         %         ylim=1-(sqrt(R^2-(xnew-Ar).^2)+B);%circle
    %         %         if ynew<ylim
    %         %             t=t-1;
    %         %         end
    %         
    %             end
    % 
    %         end
        
    end

        %immaginary check
    %         if imag(xnew)~=0 || imag(ynew)~=0
    %             1
    %              xnew
    %              ynew
    %              t
    %              xold
    %              yold
    %              xp(n,t-2)
    %              yp(n,t-2)
    %              distx
    %              disty
    %              xp(n,1)
    %              yp(n,1)
    %              pause
    %         end

        if t~=0     
            xp(n,t)=xnew;
            yp(n,t)=ynew;
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
