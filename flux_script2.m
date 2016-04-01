%computes Flux per bin as difference between particles leaving from right boundary of the bin and entering the from right boundary of the bin
%execute only after generating the volumes of the bins Vav (see scrip_stoc_diff.m)

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
                        for k=ix:P(i,it-1)-1
                            F(k,it)=F(k,it)-1;
                        end

                    elseif P(i,it-1)<ix
                        for k=P(i,it-1):ix-1
                            F(k,it)=F(k,it)+1;
                        end
                    end
                end
            end
        end
    end
end
S=zeros(T,Nx);
for ix=1:Nx
    S(:,ix)=Vav(:,ix)/(x(ix+1)-x(ix));%bin hight 
end
F=F./S';
F=F/dt;
surf(xxx,t,F(1:end-1,:)','Edgecolor','None')
 DD=-F./Cnd5;% find DD (1D diffusion) according to Fick's law (here Cnd5 is a space derivative in x) 