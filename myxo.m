%na vector with #aggregates of size from 1 to length(na) in pole A
%nb vector with #aggregates of size from 1 to length(nb) in pole B
%lambda constant release rate to incorporate in poisson distribution
%D constant diffusion step to incorporate in random walk
%dt time step
%dx space step for diffucion compartment C
%L length of diffusion compartment C

function [da,db,dc]=myxo(na,nb,lambda,D,dt,T,dx,L)
t=0:dt:T;
nt=length(t);
x=0:dx:L;
nx=length(x);
N=0;
lna=length(na);
lnb=length(nb);
nc=[];
k=0;
for i=1:lna
    for j=1:na(i)
        nc{k+j}=[0,i];
    end
    k=k+na(i);
end
for i=1:lnb
    for j=1:nb(i)
        nc{k+j}=[L,i];
    end
    k=k+nb(i);
end
da=na;
db=nb;
dc=nc;


end