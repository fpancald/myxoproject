function [xx,yy]=reflectbc2(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp)
format long
check=0;
dxx=dx/1000;
v0=[-(xnew-xold),-(ynew-yold)];
xnew0=xnew;
ynew0=ynew;
xold0=xold;
yold0=yold;

count=0;
countmax=50;
while check==0
        [xx,yy]=truncircsh2(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dxx/100);
        
%         hold on
%         plot([xold,xnew],[yold,ynew]);
%         plot([xold,xx],[yold,yy]);
%         scatter(xx,yy);

        dist=sqrt((xnew-xx)^2+(ynew-yy)^2)*damp;
        
        v=[-(xx-xold),-(yy-yold)];
        if (xx<x1+dx || xx>x2-dx) && (yy>dx || yy<w-dx)
            l=[1,0];
%             0
%             plot([xx-dxx,xx+dxx],[yy,yy]);
%             scatter(xx-dxx,yy)
%             scatter(xx+dxx,yy)
        elseif yy>=w/2 && xx>x1 && xx<=xm
            yx=(sqrt(r^2-(xx-a).^2)+b)/ynorm;
            yxdx=(sqrt(r^2-(xx-dxx-a).^2)+b)/ynorm;%circle
            l=[yxdx-yx,+dxx];
            
%             1
%             plot([xx,xx-dxx],[yx,yxdx]);
%             scatter(xx-dxx,yxdx)
        elseif yy<=w/2 && xx>x1 && xx<=xm
            yx=w-(sqrt(r^2-(xx-a).^2)+b)/ynorm;%circle
            yxdx=w-(sqrt(r^2-(xx-dxx-a).^2)+b)/ynorm;
            l=[yxdx-yx,+dxx];
%             2
%             plot([xx,xx-dxx],[yx,yxdx]);
%             scatter(xx-dxx,yxdx)
        elseif yy>=w/2 && xx>=xm && xx<x2
            yx=(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;%circle
            yxdx=(sqrt(r^2-((2*xm-xx-dxx)-a).^2)+b)/ynorm;
            l=[yxdx-yx,-dxx];
%             3
%             plot([xx,xx+dxx],[yx,yxdx]);
%             scatter(xx+dxx,yxdx)    
        elseif yy<=w/2 && xx>=xm && xx<x2
            yx=w-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;%circle
            yxdx=w-(sqrt(r^2-((2*xm-xx-dxx)-a).^2)+b)/ynorm;
            l=[yxdx-yx,-dxx];
%             4
%             plot([xx,xx+dxx],[yx,yxdx]);
%             scatter(xx+dxx,yxdx) 
        end

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
        if (x1-xx>dxx || xx-x2>dxx || yy-w>dxx ||-yy>dxx || (xx<=xm && xx>x1 && yy>w/2 && yy-(sqrt(r^2-(xx-a).^2)+b)/ynorm>dxx) || (xx<=xm && xx>x1 && yy<w/2 && -yy+w-(sqrt(r^2-(xx-a).^2)+b)/ynorm>dxx) || (xx>=xm && xx<x2 && yy>w/2 && yy-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm>dxx) || (xx>=xm && xx<x2 && yy<w/2 && -yy+w-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm>dxx))==1
            xnew=xx;
            ynew=yy;
            xold=xt;
            yold=yt;
%             v0=[-(xnew-xold),-(ynew-yold)];%need to chack
            count=count+1;
            if count>countmax
% %                 close all
                figure(11)
                hold on
                xx
                yy
                xt
                yt
                xnew0
                ynew0
                xold0
                yold0
                X1=x1:0.01:xm;
                a
                b
                r
                ynorm
        X2=xm:0.01:x2;
        Y1=(sqrt(r^2-(X1-a).^2)+b)/ynorm;
        Y2=w-(sqrt(r^2-(X1-a).^2)+b)/ynorm;
        Y3=(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
        Y4=w-(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;

        plot(X1,Y1)
        plot(X1,Y2)
        plot(X2,Y3)
        plot(X2,Y4)
        scatter(xnew,ynew,'MarkerEdgeColor','b');
        scatter(xold,yold,'MarkerEdgeColor','r');
        plot([xt,xx],[yt,yy]);
        plot([xt+500*l(1),xt-500*l(1)],[yt+500*l(2),yt-500*l(2)]);
        scatter(xx,yy,'MarkerEdgeColor','g');
        scatter(a,b);
        daspect([1 1 1])
        axis([x1 x2 0 w])
%                 [xx,yy]=reflectbc2(a,b,r,ynorm,x1,x2,xm,w,L,xnew0,ynew0,xold0,yold0,dx/10,damp);
                pause
                script_test_sing_sim
                pause

                check=1;
            end
        else
            check=1;
        end
end


end