function [xx,yy]=reflectbc(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dx,damp)
format long
check=0;
dxx=dx/1000;
v0=[-(xnew-xold),-(ynew-yold)];
while check==0
        [xx,yy]=truncircsh2(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dxx/100);
        
        hold on
        plot([xold,xnew],[yold,ynew]);
        plot([xold,xx],[yold,yy]);
        scatter(xx,yy);

        dist=sqrt((xnew-xx)^2+(ynew-yy)^2)*damp;
        
        v=[-(xx-xold),-(yy-yold)];
        if (xx<x1+dx || xx>x2-dx) && (yy>dx || yy<w-dx)
            l=[1,0];
            0
            plot([xx-dxx,xx+dxx],[yy,yy]);
            scatter(xx-dxx,yy)
            scatter(xx+dxx,yy)
        elseif yy>=w/2 && xx>x1 && xx<=xm
            yx=(sqrt(r^2-(xx-a).^2)+b)/ynorm;
            yxdx=(sqrt(r^2-(xx-dxx-a).^2)+b)/ynorm;%circle
            l=[yxdx-yx,+dxx];
            
            1
            plot([xx,xx-dxx],[yx,yxdx]);
            scatter(xx-dxx,yxdx)
        elseif yy<=w/2 && xx>x1 && xx<=xm
            yx=w-(sqrt(r^2-(xx-a).^2)+b)/ynorm;%circle
            yxdx=w-(sqrt(r^2-(xx-dxx-a).^2)+b)/ynorm;
            l=[yxdx-yx,+dxx];
            2
            plot([xx,xx-dxx],[yx,yxdx]);
            scatter(xx-dxx,yxdx)
        elseif yy>=w/2 && xx>=xm && xx<x2
            yx=(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;%circle
            yxdx=(sqrt(r^2-((2*xm-xx-dxx)-a).^2)+b)/ynorm;
            l=[yxdx-yx,-dxx];
            3
            plot([xx,xx+dxx],[yx,yxdx]);
            scatter(xx+dxx,yxdx)    
        elseif yy<=w/2 && xx>=xm && xx<x2
            yx=w-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm;%circle
            yxdx=w-(sqrt(r^2-((2*xm-xx-dxx)-a).^2)+b)/ynorm;
            l=[yxdx-yx,-dxx];
            4
            plot([xx,xx+dxx],[yx,yxdx]);
            scatter(xx+dxx,yxdx) 
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
        
        % %plot check

        X1=x1:0.001:xm;
        X2=xm:0.001:x2;
        Y1=(sqrt(r^2-(X1-a).^2)+b)/ynorm;
        Y2=w-(sqrt(r^2-(X1-a).^2)+b)/ynorm;
        Y3=(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
        Y4=w-(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;

        plot(X1,Y1)
        plot(X1,Y2)
        plot(X2,Y3)
        plot(X2,Y4)
        
        plot([x1,x1],[0,w])

        plot([x2,x2],[0,w])
        scatter(xnew,ynew);
        scatter(xold,yold,'MarkerEdgeColor','r');
        plot([xt,xx],[yt,yy]);
        plot([xt+500*l(1),xt-500*l(1)],[yt+500*l(2),yt-500*l(2)]);
        scatter(xx,yy);
        scatter(a,b);
        daspect([1 1 1])
        axis([x1 x2 0 w])
        hold off
        if (x1-xx>dxx || xx-x2>dxx || yy-w>dxx ||-yy>dxx || (xx<=xm && xx>x1 && yy>w/2 && yy-(sqrt(r^2-(xx-a).^2)+b)/ynorm>dxx) || (xx<=xm && xx>x1 && yy<w/2 && -yy+w-(sqrt(r^2-(xx-a).^2)+b)/ynorm>dxx) || (xx>=xm && xx<x2 && yy>w/2 && yy-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm>dxx) || (xx>=xm && xx<x2 && yy<w/2 && -yy+w-(sqrt(r^2-((2*xm-xx)-a).^2)+b)/ynorm>dxx))==1
            xnew=xx;
            ynew=yy;
            xold=xt;
            yold=yt;
%             v0=[-(xnew-xold),-(ynew-yold)];%need to chack
        else
            check=1;
        end
end


end