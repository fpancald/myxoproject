FF=zeros(Nx,T);
for it=2:T
    FF(:,it)=(F(:,it)-F(:,it-1));
end
% S=zeros(T,Nx);
% for ix=1:Nx
%     S(:,ix)=Vav(:,ix)/(x(ix+1)-x(ix));
% end
% FF=FF./S';
% FF=FF/dt;