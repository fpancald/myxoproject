%should load data from simulation before running this script

tic
% x=x1:(x2-x1)/Nx:x2;%close to center
x=0:(x2-x1)/Nx:L;%all domain

t=1:T;
nx=length(x)-1;
C=zeros(nx,T);

ymin=ones(T,nx)*w/2;
ymax=ones(T,nx)*w/2;

for k=1:K
    xp=XP{k};
    
    for it=1:T
        for i=1:N
            for ix=1:nx
                if xp(i,it)>=x(ix) && xp(i,it)<x(ix+1)
                    C(ix,it)=C(ix,it)+1;
                    if ymin(it,ix)>yp(i,it) 
                        ymin(it,ix)=yp(i,it);
                    end
                    if ymax(it,ix)<yp(i,it) 
                        ymax(it,ix)=yp(i,it);
                    end
                    break;%when particle is found pass to next one
                end
            end
        end
    end
    
end
C=C/K;

xx=zeros(1,nx);
V=zeros(T,nx);
Vav=zeros(T,nx);
Cnorm=zeros(nx,T);

for ix=1:nx
    xx(ix)=(x(ix)+x(ix+1))/2;
    V(:,ix)=(x(ix+1)-x(ix))*(ymax(:,ix)-ymin(:,ix));
end
for ix=1:nx
    Vav(:,ix)=(V(:,ix)+V(:,end-ix+1))/2;
    for it=1:T
        Cnorm(ix,it)=C(ix,it)/Vav(it,ix);
    end
end
der5point
cd teststat2d
mkdir debugtestalldom9
cd debugtestalldom9
save('matlab.mat','-v7.3')
toc
cd ..
cd ..