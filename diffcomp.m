function D=diffcomp(xp,T,N,Nx,L,x1,x2,dt)
xps=xp(:,2:T);
dxp=xps-xp(:,1:T-1);
% x=0:(x2-x1)/Nx:L;
x=x1:(x2-x1)/Nx:x2;
for i=1:length(x)-1
    for t=1:T-1
        dx{i,t}=[];
        for n=1:N
            if xp(n,t)>=x(i) && xp(n,t)<=x(i+1)
                dx{i,t}=[dx{i,t} dxp(n,t)];
            end
        end
    end
end
D=zeros(length(x)-1,T-1);
for i=1:length(x)-1
    for t=1:T-1
        S=var(dx{i,t});
%         S=sum(dx{i,t}.^2)/length(dx{i,t});
        D(i,t)=S/(2*dt);
    end
end
end