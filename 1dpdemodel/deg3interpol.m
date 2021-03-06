function [abcd1]=deg3interpol(x1,x2,D,Dt)
%     function [abcd1,abcd2]=deg3interpol(x1,x2,x3,D,Dt)
M1=[3*x1^2,2*x1,1,0
    3*x2^2,2*x2,1,0
    x1^3,x1^2,x1,1
    x2^3,x2^2,x2,1];
v1=[0,0,D,Dt]';
% nstep=100;
abcd1=M1\v1;
% X1=x1:(x2-x1)/nstep:x2;
% a1=abcd1(1);
% b1=abcd1(2);
% c1=abcd1(3);
% d1=abcd1(4);
% Y1=a1*X1.^3+b1*X1.^2+c1*X1+d1*ones(1,length(X1));
% Yp1=3*a1*X1.^2+2*b1*X1+c1*ones(1,length(X1));

% M2=[3*x3^2,2*x3,1,0
%     3*x2^2,2*x2,1,0
%     x3^3,x3^2,x3,1
%     x2^3,x2^2,x2,1];
% v2=[0,0,D,Dt]';
% nstep=100;
% abcd2=M2\v2;
% X2=x2+(x3-x2)/nstep:(x3-x2)/nstep:x3;
% a2=abcd2(1);
% b2=abcd2(2);
% c2=abcd2(3);
% d2=abcd2(4);
% Y2=a2*X2.^3+b2*X2.^2+c2*X2+d2*ones(1,length(X2));
% Yp2=3*a2*X2.^2+2*b2*X2+c2*ones(1,length(X2));

% X=[X1,X2];
% Y=[Y1,Y2];
% Yp=[Yp1,Yp2];
% 
% figure(1)
% plot(X,Y);
% figure(2)
% plot(X,Yp);

end
