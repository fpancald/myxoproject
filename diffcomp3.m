function [D,err]=diffcomp3(XP,T,N,Nx,L,x1,x2,dt)

% x=0:(x2-x1)/Nx:L;
x=x1:(x2-x1)/Nx:x2;

% dx=cell(length(x)-1,T-1);
for i=1:length(x)-1
    for t=1:T-1
        dx{i,t}=[];
    end
end


for t=1:T-1
%         dx{i,t}=[];
    for h=1:length(XP)
        xp=XP{h};
        for n=1:N
            for i=1:length(x)-1
                if xp(n,t)>=x(i) && xp(n,t)<=x(i+1)
                    dxp=xp(n,t+1)-xp(n,t);
                    dx{i,t}=[dx{i,t} dxp];
                    break;
                end
            end
        end
    end
end


D=zeros(length(x)-1,T-1);
err=zeros(length(x)-1,T-1);
for i=1:length(x)-1
    for t=1:T-1
        S=var(dx{i,t});
%         S=sum(dx{i,t}.^2)/length(dx{i,t});
        err(i,t)=std((dx{i,t}-mean(dx{i,t})).^2-S)/(2*dt);
        D(i,t)=S/(2*dt);
    end
end

end