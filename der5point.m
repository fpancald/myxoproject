%5-point derivative
% xxx=zeros(1,nx-1);
Cnd5=zeros(nx-1,T);
Cd5=zeros(nx-1,T);
h=x(2)-x(1);
for ix=1:nx
    
    for it=1:T
        if ix==1
            Cnd5(ix,it)=(Cnorm(ix+1,it)-Cnorm(ix,it))/(h);
            Cd5(ix,it)=(C(ix+1,it)-C(ix,it))/(h);
        elseif ix==2
            Cnd5(ix,it)=(Cnorm(ix+1,it)-Cnorm(ix-1,it))/(2*h);
            Cd5(ix,it)=(C(ix+1,it)-C(ix-1,it))/(2*h);
        elseif ix==nx
            Cnd5(ix,it)=(Cnorm(ix,it)-Cnorm(ix-1,it))/(h);
            Cd5(ix,it)=(C(ix,it)-C(ix-1,it))/(h);
        elseif ix==nx-1
            Cnd5(ix,it)=(Cnorm(ix+1,it)-Cnorm(ix-1,it))/(2*h);
            Cd5(ix,it)=(C(ix+1,it)-C(ix,it))/(2*h);
        else
            Cnd5(ix,it)=(-Cnorm(ix+2,it)+8*Cnorm(ix+1,it)-8*Cnorm(ix-1,it)+Cnorm(ix-2,it))/(12*h);
            Cd5(ix,it)=(-C(ix+2,it)+8*C(ix+1,it)-8*C(ix-1,it)+C(ix-2,it))/(12*h);
        end
    end
end