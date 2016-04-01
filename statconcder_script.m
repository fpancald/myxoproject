xd=zeros(1,Nx-1);
nCnd=zeros(Nx-1,T);
nCd=zeros(Nx-1,T);
for ix=1:Nx-1
    xd(ix)=x(ix+1);
    for it=1:T
        nCnd(ix,it)=(Cnorm(ix+1,it)-Cnorm(ix,it))/(2*(xx(ix+1)-xx(ix)));
        nCd(ix,it)=(C(ix+1,it)-C(ix,it))/(2*(xx(ix+1)-xx(ix)));
    end
end