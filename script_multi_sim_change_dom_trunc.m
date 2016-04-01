D=0.025;
x1=4.5;
x2=5.5;
xm=5;
L=10;
w=1;
T=100;
x=x1:(x2-x1)/Nx:x2;%close to center
% x=0:L/Nx:L;%all somain
t=1:T;
C=zeros(Nx,T);
ymin=ones(T,Nx)*w/2;
ymax=ones(T,Nx)*w/2;
state=1;
for k=1:K
    [xp,yp,dt]=stat_2d_diff_romr3(D,x1,x2,xm,w,N,T,L,Nx,state);

    for it=1:T
        for i=1:N
            for ix=1:Nx
                if xp(i,it)>=x(ix) && xp(i,it)<x(ix+1)
                    C(ix,it)=C(ix,it)+1;
                    if ymin(it,ix)>yp(i,it) 
                        ymin(it,ix)=yp(i,it);
                    end
                    if ymax(it,ix)<yp(i,it) 
                        ymax(it,ix)=yp(i,it);
                    end
                end
            end
        end
    end
    
end
C=C/K;

xx=zeros(1,Nx);
V=zeros(T,Nx);
Vav=zeros(T,Nx);
Cnorm=zeros(Nx,T);

for ix=1:Nx
    xx(ix)=(x(ix)+x(ix+1))/2;
    V(:,ix)=(x(ix+1)-x(ix))*(ymax(:,ix)-ymin(:,ix));
end
for ix=1:Nx
    Vav(:,ix)=(V(:,ix)+V(:,end-ix+1))/2;
    for it=1:T
        Cnorm(ix,it)=C(ix,it)/Vav(it,ix);
    end
end
der5point
