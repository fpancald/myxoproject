clear all
L=10;
rpol=1;
T=1200;
x=linspace(0,L,501);
t=linspace(0,T,5);
N0=400;
icE=85;
%u=0;
DR=1;
abcdl=deg3interpol(rpol/2,rpol,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%near left pole
abcdr=deg3interpol(L-rpol,L-rpol/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%near right pole
abcdlc=deg3interpol((L-rpol)/2,L/2,0,1);%deg3 pol (value=0 on left 1 on right der both 0)%center left 
abcdrc=deg3interpol(L/2,(L+rpol)/2,1,0);%deg3 pol (value=1 on left 0 on right der both 0)%center right
for i=1:length(t)
    for j=1:length(x)
        if  x(j)>(L-rpol)/2 && x(j)<(L+rpol)/2 %sp. diffusion and division ON near center
            Ft=(t(i))/(T);
            if x(j)<L/2
                abcd=-abcdlc*Ft;
            else
                abcd=-abcdrc*Ft;
            end
            DRtx(i,j)=(abcd(1)*x(j)^3+abcd(2)*x(j)^2+abcd(3)*x(j)+abcd(4))*DR+DR;
            DRtxDx(i,j)=(3*abcd(1)*x(j)^2+2*abcd(2)*x(j)+abcd(3))*DR;
        else

            DRtx(i,j)=DR;
            DRtxDx(i,j)=0;
        end
    end
end
figure(1)
plot(x,DRtx)
figure(2)
plot(x,DRtxDx)