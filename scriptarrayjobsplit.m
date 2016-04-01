%single simulation test script
D=0.025;
x1=4.5;
x2=5.5;
xm=5;
L=10;
w=1;
T=100;
N=10000;
Nx=100;
K=35;
% x=x1:(x2-x1)/Nx:x2;%close to center
x=0:(x2-x1)/Nx:L;%all domain

t=1:T;
nx=length(x)-1;
C=zeros(nx,T);
ymin=ones(T,nx)*w/2;
ymax=ones(T,nx)*w/2;
state=1;
for k=1:K
%     [xp,yp,dt]=stat_2d_diff_romr5(D,x1,x2,xm,w,N,T,L,Nx,state); %old version
     [xp,yp,dt]=stat_2d_diff_romr6(D,x1,x2,xm,w,N,T,L,Nx,state);% tested no receptors
%    [xp,yp,dt]=stat_2d_diff_romr7(D,x1,x2,xm,w,N,T,L,Nx,state);%1st instance receptors
    for it=1:T
        for i=1:N
            ix=locator1d(xp(i,it),0,L,(x2-x1)/Nx);
            C(ix,it)=C(ix,it)+1;
            if ymin(it,ix)>yp(i,it) 
                ymin(it,ix)=yp(i,it);
            end
            if ymax(it,ix)<yp(i,it) 
                ymax(it,ix)=yp(i,it);
            end
        end
    end
    XPs{k}=xp;
end
C=C/K;
xx=zeros(1,nx);
V=zeros(T,nx);

for ix=1:nx
    xx(ix)=(x(ix)+x(ix+1))/2;
    V(:,ix)=(x(ix+1)-x(ix))*(ymax(:,ix)-ymin(:,ix));
end

filenamecv=strcat(strcat('CV',num2str(Ntask)),'.mat');
filenamexp=strcat(strcat('XP',num2str(Ntask)),'.mat');
cd teststat2d
mkdir difftestrecsplit
cd difftestrecsplit
save(filenamecv,'C','V','-v7.3');
save(filenamexp,'XPs','-v7.3');

if Ntask==100
    flag=0;
    while flag==0
        for tag=1:99
            filename=strcat(strcat('XP',num2str(tag)),'.mat');
            if exist(filename,'file')==0
                break;
            elseif tag==99
                flag=1;
            end
        end
%         pause(1);
    end

    CC=C;
    VV=V;
    XP=XPs;
    for tag=1:99
        filename=strcat(strcat('CV',num2str(tag)),'.mat');
        load(filename);
        filename=strcat(strcat('XP',num2str(tag)),'.mat');
        load(filename);
        CC=CC+C;
        VV=VV+V;
        XP=[XP XPs];
    end
    C=CC/100;
    V=VV/100;

    Vav=zeros(T,nx);
    Cnorm=zeros(nx,T);


    for ix=1:nx
        Vav(:,ix)=(V(:,ix)+V(:,end-ix+1))/2;
        for it=1:T
            Cnorm(ix,it)=C(ix,it)/Vav(it,ix);
        end
    end
    
    cd ..
    cd ..
    der5point
    DD=diffcomp4(XP,T,N,Nx,L,x1,x2,dt);
    cd teststat2d
    cd difftestrecsplit
    
    filename=strcat(strcat('DD',num2str(Ntask)),'.mat');
    save(filename,'DD','-v7.3');
    save('data.mat','-v7.3');
end
