function [xx,yy]=truncircsh(a,b,r,x1,x2,xm,w,L,xnew,ynew,xold,yold)
tol=1e-6;
xtol=1e-3;
N=ceil(max(abs(xnew-xold),abs(ynew-yold))/xtol);
m=200;

x=xnew:(xold-xnew)/N:xold-(m-1)*(xnew-xold);
y=ynew:(yold-ynew)/N:yold-(m-1)*(ynew-yold);
% prod(y-w)
% (y-w).*y
i=2;
while (x1-x(i)>tol || x(i)-x2>tol || y(i)-w>tol ||-y(i)>tol)
i=i+1;
end
xx=x(i);
yy=y(i);

end