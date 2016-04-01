function [ix,jx] = locator2d(x,xmin,xmax,xstep,y,ymin,ymax,ystep)
%locator1d Find index location in partitioned interval
%   This function returns, based on the space coordinate x, the index
%   correspondent to the bin of the partitioned interval (Note each bin
%   contains left boundary, if x=xmax is located in separate bin)
ix=locator1d(x,xmin,xmax,xstep);
jx=locator1d(y,ymin,ymax,ystep);
end

