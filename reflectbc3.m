function [xx,yy]=reflectbc3(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp)
format long
check=0;
dxx=dx/1000;
v0=[-(xnew-xold),-(ynew-yold)];

%turn ON if control in number of reflection needed
% xnew0=xnew;
% ynew0=ynew;
% xold0=xold;
% yold0=yold;

B=w/2;
R=w/2;
Al=w/2;
Ar=L-w/2;

%turn ON if control in number of reflection needed
% count=0;
% countmax=100;
while check==0
        [xx,yy]=truncircsh3(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dxx/100);
        
%         hold on
%         plot([xold,xnew],[yold,ynew]);
%         plot([xold,xx],[yold,yy]);
%         scatter(xx,yy);

        dist=sqrt((xnew-xx)^2+(ynew-yy)^2)*damp;
        
        v=[-(xx-xold),-(yy-yold)];
        if xx<Al 
            if yy>w/2
                yx=(sqrt(R^2-(xx-Al).^2)+B);
                yxdx=(sqrt(R^2-(xx+dxx-Al).^2)+B);
                l=[yxdx-yx,-dxx];
            elseif yy<w/2
                yx=w-(sqrt(R^2-(xx-Al).^2)+B);
                yxdx=w-(sqrt(R^2-(xx+dxx-Al).^2)+B);
                l=[yxdx-yx,-dxx];
            else
                l=[1,0];
            end
            
        elseif xx>Ar
            if yy>w/2
                yx=(sqrt(R^2-(xx-Ar).^2)+B);
                yxdx=(sqrt(R^2-(xx-dxx-Ar).^2)+B);
                l=[yxdx-yx,dxx];
            elseif yy<w/2
                yx=w-(sqrt(R^2-(xx-Ar).^2)+B);
                yxdx=w-(sqrt(R^2-(xx-dxx-Ar).^2)+B);
                l=[yxdx-yx,dxx];
            else
                l=[1,0];
            end
            
        elseif (xx>=Al && xx<=x1) || (xx>=x2 && xx<=Ar)
            l=[0,1];
            
        elseif xx>x1 && xx<=xm
            if yy>w/2;
                yx=(sqrt(r^2-(xx-a).^2)+b)/ynorm;
                yxdx=(sqrt(r^2-(xx-dxx-a).^2)+b)/ynorm;%circle
                l=[yxdx-yx,+dxx];
            elseif yy<w/2 
                yx=w-(sqrt(r^2-(xx-a).^2)+b)/ynorm;%circle
                yxdx=w-(sqrt(r^2-(xx-dxx-a).^2)+b)/ynorm;
                l=[yxdx-yx,+dxx];
            else 
                l=[1,0];
            end
        elseif xx>xm && xx<x2
            if yy>w/2
                yx=(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;%circle
                yxdx=(sqrt(r^2-((2*xm-xx-dxx)-a).^2)+b)/ynorm;
                l=[yxdx-yx,-dxx];
            elseif yy<w/2 
                yx=w-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;%circle
                yxdx=w-(sqrt(r^2-((2*xm-xx-dxx)-a).^2)+b)/ynorm;
                l=[yxdx-yx,-dxx];
            else
                l=[0,1];
            end

        end
        
        %immaginary check
%         if ~isreal(l)
%             xx
%             pause
%         end

        vr=2*(v*l')/(l*l')*l-v;
        if norm(vr)~=0;
            vr=vr/norm(vr);
        else
            vr=v0;
            vr=vr/norm(vr);
        end
        xt=xx;
        yt=yy;
        xx=xx+dist*vr(1);
        yy=yy+dist*vr(2);
        
%         % %plot check
% 
%         X1=x1:0.01:xm;
%         X2=xm:0.01:x2;
%         Y1=(sqrt(r^2-(X1-a).^2)+b)/ynorm;
%         Y2=w-(sqrt(r^2-(X1-a).^2)+b)/ynorm;
%         Y3=(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
%         Y4=w-(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
% 
%         plot(X1,Y1)
%         plot(X1,Y2)
%         plot(X2,Y3)
%         plot(X2,Y4)
%         scatter(xnew,ynew);
%         scatter(xold,yold,'MarkerEdgeColor','r');
%         plot([xt,xx],[yt,yy]);
%         plot([xt+500*l(1),xt-500*l(1)],[yt+500*l(2),yt-500*l(2)]);
%         scatter(xx,yy);
%         scatter(a,b);
%         daspect([1 1 1])
%         axis([x1 x2 0 w])
%         hold off

    if xx<0 || xx>L
            check=0;
    else
        if xx<Al
             ymax=(sqrt(R^2-(xx-Al).^2)+B);
             ymin=w-(sqrt(R^2-(xx-Al).^2)+B);
         elseif xx>Ar
             ymax=(sqrt(R^2-(xx-Ar).^2)+B);
             ymin=w-(sqrt(R^2-(xx-Ar).^2)+B);
        elseif xx>x1 && xx<=xm
            ymax=(sqrt(r^2-(xx-a).^2)+b)/ynorm;
            ymin=w-(sqrt(r^2-(xx-a).^2)+b)/ynorm;
        elseif xx<x2 && xx>xm
            ymax=(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;
            ymin=w-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm; 
        else
             ymin=0;
             ymax=w;
        end
        
        if yy>ymax || yy<ymin
            check=0;
        else 
            check=1;
        end
    end

            xnew=xx;
            ynew=yy;
            xold=xt;
            yold=yt;
            
%             v0=[-(xnew-xold),-(ynew-yold)];%need to check

%control on number of reflections
%             count=count+1;
%             if count>countmax
% % %                 close all
%                 figure(11)
%                 hold on
%                 xx
%                 yy
%                 xt
%                 yt
%                 xnew0
%                 ynew0
%                 xold0
%                 yold0
%                 X1=x1:0.01:xm;
%                 a
%                 b
%                 r
%                 ynorm
%         X2=xm:0.01:x2;
%         Y1=(sqrt(r^2-(X1-a).^2)+b)/ynorm;
%         Y2=w-(sqrt(r^2-(X1-a).^2)+b)/ynorm;
%         Y3=(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
%         Y4=w-(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
% 
%         plot(X1,Y1)
%         plot(X1,Y2)
%         plot(X2,Y3)
%         plot(X2,Y4)
%         scatter(xnew,ynew,'MarkerEdgeColor','b');
%         scatter(xold,yold,'MarkerEdgeColor','r');
%         plot([xt,xx],[yt,yy]);
%         plot([xt+500*l(1),xt-500*l(1)],[yt+500*l(2),yt-500*l(2)]);
%         scatter(xx,yy,'MarkerEdgeColor','g');
%         scatter(a,b);
%         daspect([1 1 1])
%         axis([x1 x2 0 w])
% %                 [xx,yy]=reflectbc2(a,b,r,ynorm,x1,x2,xm,w,L,xnew0,ynew0,xold0,yold0,dx/10,damp);
%                 pause
%                 countmax=0;
%             end
        
end


end