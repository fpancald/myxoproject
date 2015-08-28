function [xx,yy]=truncircsh(a,b,r,ynorm,x1,x2,xm,w,L,xnew,ynew,xold,yold)
format long
tol=1e-6;
xtol=1e-3;
N=ceil(max(abs(xnew-xold),abs(ynew-yold))/xtol);
m=200;

if abs(xnew-xold)<xtol
    x=xold*ones(1,N*m+1);
else
    x=xnew:(xold-xnew)/N:xold-(m-1)*(xnew-xold);
end
if abs(ynew-yold)<xtol
    y=yold*ones(1,N*m+1);
else
    y=ynew:(yold-ynew)/N:yold-(m-1)*(ynew-yold);
end

i=2;
% %plot check
% close all
% X1=x1:0.01:xm;
% X2=xm:0.01:x2;
% Y1=(sqrt(r^2-(X1-a).^2)+b)/ynorm;
% Y2=w-(sqrt(r^2-(X1-a).^2)+b)/ynorm;
% Y3=(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
% Y4=w-(sqrt(r^2-((2*xm-X2)-a).^2)+b)/ynorm;
% figure(20)
% hold on
% plot(X1,Y1)
% plot(X1,Y2)
% plot(X2,Y3)
% plot(X2,Y4)
% scatter(xnew,ynew);
% scatter(xold,yold,'MarkerEdgeColor','r');
% plot([xold,xnew],[yold,ynew]);

while (x1-x(i)>tol || x(i)-x2>tol || y(i)-w>tol ||-y(i)>tol || (x(i)<xm && x(i)>x1 && y(i)>w/2 && y(i)-(sqrt(r^2-(x(i)-a).^2)+b)/ynorm>tol) || (x(i)<xm && x(i)>x1 && y(i)<w/2 && -y(i)+w-(sqrt(r^2-(x(i)-a).^2)+b)/ynorm>tol) || (x(i)>xm && x(i)<x2 && y(i)>w/2 && y(i)-(sqrt(r^2-((2*xm-x(i))-a).^2)+b)/ynorm>tol) || (x(i)>xm && x(i)<x2 && y(i)<w/2 && -y(i)+w-(sqrt(r^2-((2*xm-x(i))-a).^2)+b)/ynorm>tol))
%     plot([xold,x(i)],[yold,y(i)]);
    i=i+1;
end
xx=x(i);
yy=y(i);

% scatter(xx,yy);
% 
% hold off
% xx=xnew;
% yy=ynew;
end