function [Lr,Rr,rstep] = receptors(perc,N,w)
%receptors positions receptors in two matrices 
%   Matrices are filled with number of receptors per each grid block
Nr=perc*N;
dim=round(sqrt(Nr));
if mod(dim,2)~=0;
    dim=dim+1;
end
rstep=w/dim;
dimx=dim/2;
dimy=dim;
Lr=ones(dimx,dimy);
Rr=ones(dimx,dimy);
end

