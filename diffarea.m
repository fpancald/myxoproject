function [u]=diffarea()
global A0 T D
format long
A0=10;
D=0.01;
L = 10;                    %space domain size
xsteps=floor(25*L+1);       %#space points
T=200;                    %time domain size
tsteps=floor(2*T+1);        %#time points

m = 0;                      %cartesian coordinates
x = linspace(0,L,xsteps);   %space domain
t = linspace(0,T,tsteps);   %time domain

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

u = sol(:,:);

solarea = pdepe(m,@pdex1pdearea,@pdex1ic,@pdex1bc,x,t);

uarea = solarea(:,:);

figure(1)
grid off
surf(x,t,u,'EdgeColor', 'none')
title('Concentration')
xlabel('Distance x')
ylabel('Time t')

figure(2)
grid off
surf(x,t,uarea,'EdgeColor', 'none')
title('Concentration')
xlabel('Distance x')
ylabel('Time t')

end

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global D

c = 1; 
f = D*DuDx;
s=0;
end
% --------------------------------------------------------------
function [c,f,s] = pdex1pdearea(x,t,u,DuDx)
global D T

c = 1; 
f = D*DuDx;
s = -1/area(x,t)*(-D*DuDx*dxarea(x,t)+u*dtarea(x,t));

% if t<T/2
%     s = -1/area(x,t)*(-D*DuDx*dxarea(x,t)+u*dtarea(x,t));
% else
%     s = -1/area(x,t)*(-D*DuDx*dxarea(x,t));
% end
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
 if x<1
    u0=3500;
 else
     u0=0;
 end
 
% --------------------------------------------------------------
end
function [pl,ql,pr,qr] = pdex1bc(~,~,~,~,~)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end
% --------------------------------------------------------------
function A=area(x,t)
global A0 T

if x>4 && x<6 
    A=A0*(1+(x-4)*(x-6)*t/T);
else
    A=A0;
end
    
end

function dxA=dxarea(x,t)
global A0 T

if x>4 && x<6 
    dxA=A0*2*t/T*(x-5);
else
    dxA=0;
end
    
end

function dtA=dtarea(x,~)
global A0 T

if x>4 && x<6 
    dtA=A0*(x-4)*(x-6)/T;
else
    dtA=0;
end
    
end