function SmoluchowskiSolver

global V K Np chi n0

%Time parameters: NOTE: T=10^6 ==>30sec runtime
T = 5;                      %Total number of time steps.

%Particle parameters
Np = 10;                    %Max number of particles
V=[1:1:Np]';                %Volumes/sizes of particles  1 to 10 cm^3
a = ((3/(4*pi)).*V).^(1/3); %Radius corresponding to Volume

%Initial Conditions
%n0 = zeros(Np+1,1);
n0 = zeros(Np,1);
n0(1) = 10;
n0(2) = 0;
%n0(Np+1) = sum(V.*n0(1:Np));


%Compute Coagulation Kernel
for (p=1:Np)
    for(q=1:Np)
        K(p,q) = 1;%     pi*((a(p)+a(q))^2);%   V(p)*V(q);%     
    end
end

%Solve Smoluchowski Equation
[t,n] = ode45(@SmoluchowskiEqns, [0 T], n0);

%Check conservation
chi = zeros(length(t),1);
for (i=1:length(t))
    chi(i) = V(1:Np)'*n(i,1:Np)';
end

%Plot and Save resutls
h= figure;
hold on
grid on
plot(t,n,'LineWidth',2)
plot(t,chi,'--', 'LineWidth',2)
xlabel('Time t','FontSize', 16);
ylabel('Number of Particles for each Volume Size','FontSize', 16);
titletext = ['No Max Volume'];
title(titletext,'FontSize', 20);
axis([min(t), max(t), 0, max(n0)+0.5])
picfilename = ['ODE45_SmoluchowksiNoMax_Neqn' num2str(Np)]
saveas(h,picfilename,'pdf');
%close all
 
end

function dn = SmoluchowskiEqns(t,n)

global V K Np chi n0

%n(Np+1) = sum(V.*n(1:Np));

%Define Smoluchowski ODEs
%dn = zeros(Np+1,1);
dn = zeros(Np,1);
for (p=1:Np)
	nIntoP = 0;
    nOutofP = 0;
    for (q=1:Np)
    	if (q<=p-1)
        	nIntoP = nIntoP + 0.5*K(p-q,q)*n(p-q)*n(q);
        end
        nOutofP = nOutofP + K(p,q)*n(p)*n(q);
    end
    dn(p) = nIntoP - nOutofP;
end

%Enforce conservation
0*dn(1) == -V(2:Np)'*dn(2:Np)/V(1);

%0*d(Np+1) == V(1:Np)'*n(1:Np) - n0(Np+1);

end