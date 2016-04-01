function [xx,yy]=truncircsh3(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold,dxx)
format long
% tol=dxx;
xtol=dxx*10;
N=ceil(max(abs(xnew-xold),abs(ynew-yold))/xtol);
m=200;

B=w/2;
R=w/2;
Al=w/2;
Ar=L-w/2;

if abs(xnew-xold)<xtol
    x=xold*ones(1,N*m+2);
else
%     x=xnew:(xold-xnew)/N:xold-(m-1)*(xnew-xold);%from new to old
    x=xold:-(xold-xnew)/N:xnew;
end
if abs(ynew-yold)<xtol
    y=yold*ones(1,N*m+2);
else
%     y=ynew:(yold-ynew)/N:yold-(m-1)*(ynew-yold);%from new to old
    y=yold:-(yold-ynew)/N:ynew;
end


i=1;


% %plot check
% close all
% 
% X1=x1:0.01:xm;
% X2=xm:0.01:x2;
% Y1=(sqrt(r^2-(X1-a).^2)+b)/ynorm;
% Y2=w-(sqrt(r^2-(X1-a).^2)+b)/ynorm;
% Y3=(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
% Y4=w-(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
% figure(20)
% hold on
% plot([x(1),x(end)],[y(1),y(end)])
% plot(X1,Y1)
% plot(X1,Y2)
% plot(X2,Y3)
% plot(X2,Y4)
% scatter(xnew,ynew);
% scatter(xold,yold,'MarkerEdgeColor','r');
% plot([xold,xnew],[yold,ynew]);

%from new to old
% while (x1-x(i)>tol || x(i)-x2>tol || y(i)-w>tol ||-y(i)>tol || (x(i)<=xm && x(i)>x1 && y(i)>w/2 && y(i)-(sqrt(r^2-(x(i)-a).^2)+b)/ynorm>tol) || (x(i)<=xm && x(i)>x1 && y(i)<w/2 && -y(i)+w-(sqrt(r^2-(x(i)-a).^2)+b)/ynorm>tol) || (x(i)>=xm && x(i)<x2 && y(i)>w/2 && y(i)-(sqrt(r^2-((2*xm-x(i))-a).^2)+b)/ynorm>tol) || (x(i)>=xm && x(i)<x2 && y(i)<w/2 && -y(i)+w-(sqrt(r^2-((2*xm-x(i))-a).^2)+b)/ynorm>tol))==1
% %     plot([xold,x(i)],[yold,y(i)]);
%     i=i+1;
% end
% xx=x(i);
% yy=y(i);

%from old to new
check=1;
while check==1
    if x(i)<0 || x(i)>L
        check=0;
    else
        if x(i)<Al
             ymax=(sqrt(R^2-(x(i)-Al).^2)+B);
             ymin=w-(sqrt(R^2-(x(i)-Al).^2)+B);
         elseif x(i)>Ar
             ymax=(sqrt(R^2-(x(i)-Ar).^2)+B);
             ymin=w-(sqrt(R^2-(x(i)-Ar).^2)+B);
        elseif x(i)>x1 && x(i)<=xm
            ymax=(sqrt(r^2-(x(i)-a).^2)+b)/ynorm;
            ymin=w-(sqrt(r^2-(x(i)-a).^2)+b)/ynorm;
        elseif x(i)<x2 && x(i)>xm
            ymax=(sqrt(r^2-((2*xm-x(i))-a).^2)+b)/ynorm;
            ymin=w-(sqrt(r^2-((2*xm-x(i))-a).^2)+b)/ynorm; 
        else
             ymin=0;
             ymax=w;
        end
        
        if y(i)>ymax || y(i)<ymin
            check=0;
        else 
            check=1;
        end
    end
        
        
%     plot([xold,x(i)],[yold,y(i)]);
if check==1
    i=i+1;
end
%     if x(i-1)==xm && (y(i-1)==(sqrt(r^2-(x(i-1)-a).^2)+b)/ynorm || y(i-1)==w-(sqrt(r^2-(x(i-1)-a).^2)+b)/ynorm)
%         
%         break;
%     end
   
    if i>min(length(x),length(y))    
        break;
    end
%     i
end
% i
% x(i)
% y(i)
xx=x(i-1);
yy=y(i-1);

%check on domain constraint
% if xx>L
%     xx
%     x(i-2)
%     pause
% end
% scatter(xx,yy);
% 
% hold off

end