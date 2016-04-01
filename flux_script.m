F=zeros(Nx,T);
P=zeros(N,T);
tol=1e-6;
for it=1:T
    for i=1:N
        for ix=1:Nx
            if xp(i,it)>=x(ix)-tol && xp(i,it)<x(ix+1)+tol
                P(i,it)=ix;
                if it>1
                    if P(i,it-1)>ix
                        F(ix,it)=F(ix,it)-1;

                    elseif P(i,it-1)<ix
                        F(P(i,it-1),it)=F(P(i,it-1),it)+1;

                    end
                end
            end
        end
    end
end
S=zeros(T,Nx);
for ix=1:Nx
    S(:,ix)=Vav(:,ix)/(x(ix+1)-x(ix));
end
F=F./S';
% DD=-F./Cnd;