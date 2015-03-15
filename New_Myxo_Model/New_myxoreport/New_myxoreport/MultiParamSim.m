global iter
niter=1; %number of simulations
for iter=1:niter
%Run parameter Sweep for RomR Diffusion

close all;
clear
%default parameter values
%load the fmax equatiom

global alpha DD Dd p0 DiffusionLimitingBigZone DiffusionLimitingNarZone sigma2 fmax CarryingCapicity a CellDiv L T N1 N2 n3 sigma1 r beta Tbd nmr
N1=1100;% number of receptors pole1
N2=600;% number of receptors pole2
n3=400;% number of receptors division site
nmr=4.5;% average maximum number of binding molecules per receptor
sigma1=0.01; % basic association rate
r=0.1;%portion of cell corresponding to pole
L = 10;% cell length
Tbd=60*10;%time in seconds before division starts 
T=60*20; %time in seconds after division starts
beta=1;% carrying capacity speed parameter
%load the fmax equatiom
%SepturmArea

%list of times with gamma distribution for reversal
a= gamrnd(1.4, 78.2, 10, 1); %reversal at a random time (Gamma dist.)

% a=60*ones(10,1); %reversal at a fixed time

for l=2:length(a)%time sequence for reversal
        a(l)=a(l)+a(l-1);
end
a=a*10; % to make time in frames rather than seconds, multiply by 10 for 10 seconds/frame

%load('fittedmodelSpline.mat');% load a carrying capacity function
%fmax = fittedmodel;
alpha = 10; %evolution parameter for division site
%Alpha_values = [10 5 1 1/5 1/10];
%alpha = Alpha_values(i)
sigma2 = 0.01; %dissociation rate
DD_values = [2.4 ];% 1.2 0.6 ];%5; %diffusion rate of free RomR
Dd_values = [1e-4 ];% 1e-4 1e-4]; %diffusion of bounded RomR
p0 = 0;%production rate
DiffusionLimitingBigZone = 0; % this limits the free RomR diffusion at the midpoint ON/OFF (1/0)
DiffusionLimitingNarZone = 0;
CellDiv=1; %cell division ON/OFF (1/0)
CarryingCapicity = 1; % carrying capacity for bound RomR ON/OFF (1/0)

for i= 1:length(DD_values)
DD = DD_values(i);
Dd = Dd_values(i);
%[RomRavgfinal u1, u2] = mixoproj4(); % non-limiting polar receptor
[RomRavgfinal, u1, u2] = mixoproj5cc(); % limiting polar receptor case
end

end