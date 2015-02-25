N=2;
dt=0.1;
% mx=[];
% for it=1:s(1)
%     mx=[mx max(D{it})];
% end
% mmx=max(mx);
% for it=1:s(1)
%     plot(X{it},D{it})
%     ylim([0,0.1*mmx])
%     pause(0.01)
% end

sm=[];
mx=[];
for it=1:s(1)
    sm=[sm sum(D{it})];
    mx=[mx max(D{it}/sm(end))];
end
% mmx=min(mx);
mmx=max(mx);
for it=1:s(1)
%     plot(X{it},D{it}/sm(it))
    plot(X{it},D{it}/sm(it)*100)
%     ylim([0,mx(end)])
%     ylim([0,mmx])
%     ylim([0,mx(end)*100])
    ylim([0,mmx*100])
    title(num2str(dt*(it-1)));
    if (mod(it,N)==1)
        picfilename=sprintf('%s%05d%s','DensityPics/pic_',it,'.png');
        saveas(gcf,picfilename);
    end

%     pause(0.01)
end

% mx=[];
% for it=1:s(1)
%     mx=[mx max(D{it})];
% end
% mmx=min(mx);
% for it=1:s(1)
%     plot(X{it},D{it}/mx(it))
%     ylim([0,1])
%     pause(0.01)
% end