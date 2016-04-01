% x=0:0.1:10;
% y=x.^2;
% z=sqrt(x);
% fig1=figure(1);
% plot(x,y);
% fig2=figure(2);
% plot(x,z);
% 
% print('-f1','plot1','-dpng');
% print('-f2','plot2','-dpng');
% 
% saveas(fig1,'plot1.fig');
% saveas(fig2,'plot2.fig');
cd teststat2d
cd 20160331T000000020
load('data.mat','C','Cnorm','Cd5','Cnd5','t','T','xx')
cd ..
cd ..
plotstat2d
cd teststat2d
cd 20160331T000000020
fig1=figure(1);
fig2=figure(2);
fig8=figure(8);
fig9=figure(9);
fig11=figure(11);
fig12=figure(12);
saveas(fig1,'plot1.fig');
saveas(fig2,'plot2.fig');
saveas(fig8,'plot8.fig');
saveas(fig9,'plot9.fig');
saveas(fig11,'plot11.fig');
saveas(fig12,'plot12.fig');

