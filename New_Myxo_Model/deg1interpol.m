function [ab1]=deg1interpol(x1,x2,N,N0)
M1=[x1,1
    x2,1];
v1=[N,N0]';
% nstep=100;
ab1=M1\v1;
% X1=x1:(x2-x1)/nstep:x2;
% a1=ab1(1);
% b1=ab1(2);
% Y1=a1*X1+b1*ones(1,length(X1));

% figure(1)
% plot(X1,Y1);


end