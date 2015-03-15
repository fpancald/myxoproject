function [s1profile] = bindingRateProfile(alpha)
%binding rate
global a CellDiv L T N1 N2 n3 sigma1 r Tbd tpts xpts
x = linspace(0,L,xpts);
Allt = linspace(0,T+Tbd,tpts);
s1profile = ones(length(x), length(Allt));
s1=sigma1;
s1profile=s1*s1profile;
for jj = 1:length(Allt)
 t = Allt(jj);   
for ii = 1:length(x) 
for irev=1:length(a)
    if t<a(irev)
        if mod(irev,2)==1
            n1=N1;
            n2=N2;
        else
            n1=N2;
            n2=N1;
        end
        break;
    end
end
if (t>Tbd)
    td=t-Tbd;
else
    td=0;
end
if (x(ii)<r*L/2)
    s1profile(ii,jj)=s1*n1;
elseif (x(ii)>=r*L/2 && x(ii)<r*L)
    s1profile(ii,jj)=((2*(1-n1))*x(ii)/(r*L)+2*n1-1)*s1;
elseif (x(ii)>L/2*(1-r) && x(ii)<=L/2)
    %s1profile(ii)=(-2/(L*r)*(1-n1)*td/T*x(ii)-(n1-1)/r*td/T+1+(n1-1)*td/T)*s1;
    s1profile(ii,jj)=(-2/(L*r)*(1-n3)*(td/T*CellDiv)^(1/alpha)*x(ii)-(n3-1)/r*(td/T*CellDiv)^(1/alpha)+1+(n3-1)*(td/T*CellDiv)^(1/alpha))*s1;
elseif (x(ii)>L/2 && x(ii)<=L/2*(1+r))
    %s1profile(ii)=(+2/(L*r)*(1-n2)*td/T*x(ii)+(n2-1)/r*td/T+1+(n2-1)*td/T)*s1;
    s1profile(ii,jj)=(+2/(L*r)*(1-n3)*(td/T*CellDiv)^(1/alpha)*x(ii)+(n3-1)/r*(td/T*CellDiv)^(1/alpha)+1+(n3-1)*(td/T*CellDiv)^(1/alpha))*s1;   
elseif (x(ii)>L*(1-r) && x(ii)<L*(1-r/2))
    s1profile(ii,jj)=s1*(2/(r*L)*(n2-1)*x(ii)+(2*n2-1-2/r*(n2-1)));
elseif(x(ii)>=L*(1-r/2))
    s1profile(ii,jj)=n2*s1;
end
end
end