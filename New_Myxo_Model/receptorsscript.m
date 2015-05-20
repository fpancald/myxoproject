clear all
L=10;
rpol=1;
T=1200;
x=linspace(0,L,501);
t=linspace(0,T,5);
N0=400;
icE=85;
%u=0;

abcdl=deg3interpol(rpol/2,rpol,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%near left pole
abcdr=deg3interpol(L-rpol,L-rpol/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%near right pole
abcdlc=deg3interpol((L-rpol)/2,L/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%center left 
abcdrc=deg3interpol(L/2,(L+rpol)/2,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%center right
for i=1:length(t)
    for j=1:length(x)
        u(i,j)=2*icE*sin(x(j)*i);
        if x(j)<rpol/2 ||x(j)>L-rpol/2 %poles
            N(i,j)=N0*(1+u(i,j)/(4*icE));%linear interaction

        elseif x(j)>=rpol/2 && x(j)<rpol %near left pole
            abcd=abcdl*N0*(1+u(i,j)/(4*icE));
            N(i,j)=abcd(1)*x(j)^3+abcd(2)*x(j)^2+abcd(3)*x(j)+abcd(4);%deg 3 receptor spline

        elseif x(j)>L-rpol && x(j)<=L-rpol/2%near right pole
            abcd=abcdr*N0*(1+u(i,j)/(4*icE));
            N(i,j)=abcd(1)*x(j)^3+abcd(2)*x(j)^2+abcd(3)*x(j)+abcd(4);

        elseif x(j)>(L-rpol)/2 && x(j)<(L+rpol)/2 %&& Div==1 && t>divT%division ON,dvision started,and near center
            Ft=(t(i))/(T);
            if x(j)<L/2%left of center

                abcd=abcdlc*N0/2*(1+u(i,j)/(4*icE))*Ft;
                N(i,j)=abcd(1)*x(j)^3+abcd(2)*x(j)^2+abcd(3)*x(j)+abcd(4);

            else%right of center
                abcd=abcdrc*N0/2*(1+u(i,j)/(4*icE))*Ft;
                N(i,j)=abcd(1)*x(j)^3+abcd(2)*x(j)^2+abcd(3)*x(j)+abcd(4);
            end

        else%not close to center,not poles or close, division OFF or division not started
            N(i,j)=0;
        end
    end
end
plot(x,N)
