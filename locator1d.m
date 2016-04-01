function [ix] = locator1d(x,xmin,xmax,xstep)
%locator1d Find index location in partitioned interval
%   This function returns, based on the space coordinate x, the index
%   correspondent to the bin of the partitioned interval (Note each bin
%   contains left boundary, if x=xmax is located in last bin)
n=length(x);
ix=zeros(1,n);
for i=1:n
    if x(i)>=xmin && x(i)<xmax
        ix(i)=floor((x(i)-xmin)/xstep)+1;
    elseif x(i)==xmax
        ix(i)=floor((x(i)-xmin)/xstep);
    end
    
end
end

