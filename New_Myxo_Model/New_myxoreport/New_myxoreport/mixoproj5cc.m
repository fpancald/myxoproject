function [RomRavgfinal, u1, u2] = mixoproj5cc()
format long
m = 0;
global alpha DD Dd p0 DiffusionLimitingBigZone DiffusionLimitingNarZone sigma2 iter L T Tbd nmr tpts xpts 
xpts=251;
tpts=1801;
x = linspace(0,L,xpts);
t = linspace(0,T+Tbd,tpts); %T+Tbd = Total time in seconds

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
AssymetryRatio = sum((u1(end, 127:151)+u2(end, 127:151))) / sum(u1(end, 101:125)+u2(end, 101:125));

%get Time lapse of sigma1 as s1profile
s1profile = bindingRateProfile(alpha);

%time lapse
mixoplot=figure('visible','off');
myColors = colormap(jet(length(t)));
for j = 1:250:length(t)
  %color timeplot
  thisColor = myColors(j, :);
%bw timeplots
  %maxWhite = .9;
  %thisColor = [max((maxWhite-maxWhite*(j/2000)), 0) max((maxWhite-maxWhite*(j/2000)), 0) max((maxWhite-maxWhite*(j/2000)), 0)];

    subplot(413)    
    plot(x,(u1(j,:)+u2(j,:))/(u1(j,1)+u2(j,1)), '-', 'color', thisColor)
    ylabel({'Norm Total'; 'RomR(x)'}, 'Fontsize', 16)
    hold on
    
    subplot(412)
    plot(x,u1(j,:),'-', 'color', thisColor)
    xlabel('Distance x')
    ylabel({'Free' ;'RomR(x)'}, 'Fontsize', 16)
    hold on
    
    subplot(411)
    title({['Time t =', num2str(t(j))];['Assymetry = ' , num2str(AssymetryRatio)]}, 'Fontsize', 20)
    plot(x,u2(j,:),'-', 'color', thisColor)
    ylabel({'Bound'; 'RomR(x)'}, 'Fontsize', 16)
    hold on
    
    subplot(414)
    plot(x,s1profile(:,j),'-', 'color', thisColor)
    ylabel('Sigma1(x)', 'Fontsize', 16)
    xlabel({['Alpha = ', num2str(alpha), ';  DD =  ', num2str(DD), '; Dd = ', num2str(Dd),...
        ' ; S2 = ' , num2str(sigma2)]; [' p0 = ' , num2str(p0), '; DiffZone2 = ', num2str(DiffusionLimitingBigZone), ...
        '; DiffZone1 = ', num2str(DiffusionLimitingNarZone)] }, 'Fontsize', 20)
    hold on

end

MixoPlotFileName = ['Images/plotFiles_mixoPlot_', num2str(iter,'%05d'),'.png'];
saveas(mixoplot, MixoPlotFileName);
for j = 1:25:length(t)
    
    mixoplot=figure('visible','off');    
    plot(x,(u1(j,:)+u2(j,:))/(u1(j,1)+u2(j,1)), '-', 'color', thisColor)
    title(['Time t =', num2str(t(j))]);
    xlabel('Distance x')
    ylabel({'Norm Total'; 'RomR(x)'}, 'Fontsize', 16)
    MixoPlotFileName = ['Images/plotFiles_norm_mixoPlot_', num2str(iter,'%05d'),'time=',num2str(t(j)),'.png'];
    saveas(mixoplot, MixoPlotFileName);

    mixoplot=figure('visible','off'); 
    plot(x,u1(j,:),'-', 'color', thisColor)
    title(['Time t =', num2str(t(j))]);
    xlabel('Distance x')
    ylabel({'Free' ;'RomR(x)'}, 'Fontsize', 16)
    ylim([0 max(max(u1))]);
    MixoPlotFileName = ['Images/plotFiles_mixoPlot_', num2str(iter,'%05d'),'time=',num2str(t(j)),'.png'];
    saveas(mixoplot, MixoPlotFileName);

    mixoplot=figure('visible','off'); 
    plot(x,u2(j,:),'-', 'color', thisColor)
    title(['Time t =', num2str(t(j))]);
    xlabel('Distance x')
    ylabel({'Bound'; 'RomR(x)'}, 'Fontsize', 16)
    ylim([0 max(max(u2))]);
    MixoPlotFileName = ['Images/plotFiles_mixoPlot_', num2str(iter,'%05d'),'time=',num2str(t(j)),'.png'];
    saveas(mixoplot, MixoPlotFileName);

    mixoplot=figure('visible','off'); 
    plot(x,s1profile(:,j),'-', 'color', thisColor)
    title(['Time t =', num2str(t(j))]);
    ylabel('Sigma1(x)', 'Fontsize', 16)
    xlabel({['Alpha = ', num2str(alpha), ';  DD =  ', num2str(DD), '; Dd = ', num2str(Dd),...
        ' ; S2 = ' , num2str(sigma2)]; [' p0 = ' , num2str(p0), '; DiffZone2 = ', num2str(DiffusionLimitingBigZone), ...
        '; DiffZone1 = ', num2str(DiffusionLimitingNarZone)] }, 'Fontsize', 20)
    ylim([0 max(max(s1profile))]);
    MixoPlotFileName = ['Images/plotFiles_mixoPlot_', num2str(iter,'%05d'),'time=',num2str(t(j)),'.png'];
    saveas(mixoplot, MixoPlotFileName);
    close all
end
 RomR=u1(1,:)+u2(1,:);
 for j = 2:length(t)
 RomR=RomR+u1(j,:)+u2(j,:);
 end

RomR=RomR/length(t);
RomRavgfinal = RomR;
end

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global alpha DD Dd p0 DiffusionLimitingBigZone DiffusionLimitingNarZone sigma2 fmax CarryingCapicity a CellDiv L T N1 N2 n3 r sigma1 beta Tbd nmr

for irev=1:length(a)
    if t<a(irev)
        if mod(irev,2)==1
            n1=N1;
            n2=N2;
            %fm = fmax(x);
        else
            n1=N2;
            n2=N1;
            %fm = fmax(10-x);
        end
        break;
    end
end
myDD = DD;
s1=sigma1;
p1 = 0; % production of unbound
p2 = 0; % proeuction of bound
if (x<r*L/2)    
    % carrying capacity 
    nr=n1;
%     if CarryingCapicity
%       s1 = s1*(1 - u(2)/fm)^(1/beta)*nr;
%     else
      s1=s1*nr;    
%     end
      p1 = p0*nr;
elseif (x>=r*L/2 && x<r*L)
    % carrying capacity
    nr=((2*(1-n1))*x/(r*L)+2*n1-1);
%     if CarryingCapicity
%       s1=nr*s1*(1 - u(2)/fm)^(1/beta);
%     else
        s1=nr*s1;
%     end
    %p1 = p0*n1;   
elseif (x>L*(1-r) && x<L*(1-r/2))
    %s1=s1*(2/(r*L)*(n2-1)*x+(2*n2-1-2/r*(n2-1)));
    % carrying capacity
    nr=(2/(r*L)*(n2-1)*x+(2*n2-1-2/r*(n2-1)));
%     if CarryingCapicity
%       s1=s1*(1 - u(2)/fm)^(1/beta)*nr;
%     else
     s1=s1*nr;        
%     end
    %p1 = p0*n2;
elseif(x>=L*(1-r/2))
    % carrying capacity
    nr=n2;
%     if CarryingCapicity
%      s1= s1*(1 - u(2)/fm)^(1/beta)*nr;
%     else
     s1=nr*s1;        
%     end
     p1 = p0*nr;
else
    if(t<Tbd)
        nr=1;
        s1=nr*s1;
        myDD=DD;
    else
        td=t-Tbd;
        if (x>L/2*(1-r) && x<=L/2*(1-r/2))
            %s1=(-2/(L*r)*(1-n1)*td/T*x-(n1-1)/r*td/T+1+(n1-1)*td/T)*s1;
            nr=(-2/(L*r)*(1-n3)*(td/T*CellDiv)^(1/alpha)*x-(n3-1)/r*(td/T*CellDiv)^(1/alpha)+1+(n3-1)*(td/T*CellDiv)^(1/alpha));
            s1=nr*s1;
            %p1 = 100*p0*n2;
             if DiffusionLimitingBigZone
                myDD=DD*(1-td/T*CellDiv)^(alpha);
               % Dd=Dd*(1-(td/T*CellDiv)^(1/alpha));
             end
        % Added A narrower diffusion limiting zone: only half of the binding zone   
        elseif (x>L/2*(1-r/2) && x<=L/2)
            %s1=(-2/(L*r)*(1-n1)*td/T*CellDiv*x-(n1-1)/r*td/T*CellDiv+1+(n1-1)*td/T*CellDiv)*s1;
            nr=(-2/(L*r)*(1-n3)*(td/T*CellDiv)^(1/alpha)*x-(n3-1)/r*(td/T*CellDiv)^(1/alpha)+1+(n3-1)*(td/T*CellDiv)^(1/alpha));
            s1=nr*s1;
            %p1 = 100*p0*n2;
            if DiffusionLimitingNarZone
               myDD=DD*(1-td/T*CellDiv)^(alpha);
              % Dd=Dd*(1-(td/T*CellDiv)^(1/alpha));
            end    
        elseif (x>L/2 && x<=L/2*(1+r/2))
            %p1 = -100*p0*n2;
            %s1=(+2/(L*r)*(1-n2)*td/T*CellDiv*x+(n2-1)/r*td/T*CellDiv+1+(n2-1)*td/T*CellDiv)*s1;
            nr=(+2/(L*r)*(1-n3)*(td/T*CellDiv)^(1/alpha)*x+(n3-1)/r*(td/T*CellDiv)^(1/alpha)+1+(n3-1)*(td/T*CellDiv)^(1/alpha));
            s1=nr*s1;
            if DiffusionLimitingNarZone
               myDD=DD*(1-td/T*CellDiv)^(alpha);
              % Dd=Dd*(1-(td/T*CellDiv)^(1/alpha));
            end
        % End of Diffusion limiting Zone    
        elseif (x>L/2*(1+r/2) && x<=L/2*(1+r))
            %p1 = -100*p0*n2;
            %s1=(+2/(L*r)*(1-n2)*td/T*CellDiv*x+(n2-1)/r*td/T*CellDiv+1+(n2-1)*td/T*CellDiv)*s1;
            nr=(+2/(L*r)*(1-n3)*(td/T*CellDiv)^(1/alpha)*x+(n3-1)/r*(td/T*CellDiv)^(1/alpha)+1+(n3-1)*(td/T*CellDiv)^(1/alpha));
            s1=nr*s1;
             if DiffusionLimitingBigZone
                myDD=DD*(1-td/T*CellDiv)^(alpha);
               % Dd=Dd*(1-(td/T*CellDiv)^(1/alpha));
             end
        else
           nr=1;
            s1=nr*s1;
            myDD=DD;
        end
    end
end
    if CarryingCapicity
        fm=nmr*nr;
      s1 = s1*(1 - u(2)/fm)^(1/beta);   
    end
%p2 = -p1;
s2= sigma2;
c = [1;1]; 
f = [myDD+Dd;Dd].*DuDx; 

% adding the terms for the constricting area
gsig = .1;
g = @(xx)(exp(-(xx-(L/2)).^2/(2*gsig^2)));
amax = @(tt)  4.9*(tt/T*CellDiv)^(1/alpha);
A0 = 5;
A = A0 - amax(t)*g(x);
A_t = -g(x)*(t / T*CellDiv)^(1/alpha) / (alpha * T);
A_x = (x - L/2)*amax(t)*g(x)/(gsig^2)*CellDiv;

s = [(-u(1)*A_t - DD*DuDx(1)*A_x)/A -s1*u(1)+s2*u(2)+p1  ;  (-u(2)*A_t - Dd*DuDx(2)*A_x)/A +s1*u(1)-s2*u(2)+p2];
s = s(:,2);
%s = [-s1*u(1)+s2*u(2)+p1;+s1*u(1)-s2*u(2)+p2];
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
global L r
f=r;

if (x<f*L/2)
    u0=[0;3000];
%     u0=[0;2300];
elseif (x<f*L)
    u0=[0;3000];
%     u0=[0;2300];
elseif(x>f*L && x<L*(1-f))
    u0=[0;50];
elseif(x>L*(1-f))
    u0=[0;1600];
%     u0=[0;2300];
else
    u0=[0;1600];
%     u0=[0;2300];
end

end
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = [0;0];
ql = [1;1];
pr = [0;0];
qr = [1;1];
end
