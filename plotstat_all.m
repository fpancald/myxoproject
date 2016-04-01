figure(1)
surf(xx,t,Cnorm','Edgecolor','None')
figure(2)
surf(xx,t,C','Edgecolor','None')

% for it=1:T
%     pause(0.1)
%     figure(3)
%     plot(xx,V(it,:))
%     axis([x1 x2 0 (x2-x1)/Nx*w]);
% end
% for it=1:T
%     pause(0.1)
%     figure(4)
%     plot(xx,Vav(it,:))
%     axis([x1 x2 0 (x2-x1)/Nx*w]);
% end
% for it=1:T
%     pause(0.1)
%     figure(5)
%     plot(xx,Cnorm(:,it));
%     axis([x1 x2 0 0.05*T*N]);
% end
% figure(6)
% plot(xx,Cnorm(:,1));
% axis([x1 x2 0 N]);
% figure(7)
% plot(xx,Cnorm(:,end));
% axis([x1 x2 0 N]);

figure(8)
surf(xx,t,Cnd5','Edgecolor','None')
figure(9)
surf(xx,t,Cd5','Edgecolor','None')
figure(12)
kk=0;
for k=0:20:T
    kk=kk+1;
    subplot(3,2,kk)
    if k~=T
        plot(xx,Cnd5(:,k+1)')
    else
        plot(xx,Cnd5(:,T)')
    end
end
figure(11)
kk=0;
for k=0:20:T
    kk=kk+1;
    subplot(3,2,kk)
    if k~=T
        plot(xx,Cnorm(:,k+1)')
    else
        plot(xx,Cnorm(:,T)')
    end
end