mx=[];
for it=1:s(1)
    mx=[mx max(D{it})];
end
mmx=max(mx);
for it=1:s(1)
    plot(X{it},D{it})
    ylim([0,0.1*mmx])
    pause(0.01)
end

% sm=[];
% mx=[];
% for it=1:s(1)
%     sm=[sm sum(D{it})];
%     mx=[mx max(D{it}/sm(end))];
% end
% % mmx=min(mx);
% mmx=max(mx);
% for it=1:s(1)
%     plot(X{it},D{it}/sm(it))
% %     ylim([0,mx(end)])
%     ylim([0,mmx])
%     pause(0.01)
% end

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