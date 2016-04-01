%single simulation test script
clear all
tic
D=0.025;
x1=4.5;
x2=5.5;
xm=5;
L=10;
w=1;
T=100;
N=10000;
Nx=100;
K=10008;
% x=x1:(x2-x1)/Nx:x2;%close to center
x=0:(x2-x1)/Nx:L;%all domain

t=1:T;
nx=length(x)-1;
C=zeros(nx,T);
Ymin=zeros(T,nx);
Ymax=zeros(T,nx);


state=1;
matlabpool 12
parfor k=1:K
%     xx=x1:(x2-x1)/Nx:x2;%close to center
    xx=0:(x2-x1)/Nx:L;%all domain
    
    ymin{k}=ones(T,nx)*w/2;
    ymax{k}=ones(T,nx)*w/2;
    CC{k}=zeros(nx,T);
    [xp,yp,dt]=stat_2d_diff_romr5(D,x1,x2,xm,w,N,T,L,Nx,state);

    for it=1:T
        for i=1:N
            for ix=1:length(xx)
                if xp(i,it)>=xx(ix) && xp(i,it)<xx(ix+1)
                    CC{k}(ix,it)=CC{k}(ix,it)+1;
                    if ymin{k}(it,ix)>yp(i,it) 
                        ymin{k}(it,ix)=yp(i,it);
                    end
                    if ymax{k}(it,ix)<yp(i,it) 
                        ymax{k}(it,ix)=yp(i,it);
                    end
                end
            end
        end
    end
%     C=C+CC{k};
end
matlabpool close


for k=1:K
    C=C+CC{k};
    Ymax=Ymax+ymax{k};
    Ymin=Ymin+ymin{k};
end
Ymax=Ymax/K;
Ymin=Ymin/K;
C=C/K;

xx=zeros(1,nx);
V=zeros(T,nx);
Vav=zeros(T,nx);
Cnorm=zeros(nx,T);

for ix=1:nx
    xx(ix)=(x(ix)+x(ix+1))/2;
    V(:,ix)=(x(ix+1)-x(ix))*(Ymax(:,ix)-Ymin(:,ix));
end
for ix=1:nx
    Vav(:,ix)=(V(:,ix)+V(:,end-ix+1))/2;
    for it=1:T
        Cnorm(ix,it)=C(ix,it)/Vav(it,ix);
    end
end
der5point
cd teststat2d
mkdir debugtestalldompar2
cd debugtestalldompar2
save('matlab.mat','-v7.3')
toc
cd ..
cd ..
