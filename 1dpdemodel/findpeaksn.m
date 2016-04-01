test=u5+u6;

% lflag(1)=0;
% rflag(1)=0;
% for it=10:length(t)
%     if test(it,1)>test(it-1,1)
%         lflag(it)=1;
%     elseif test(it,1)<test(it-1,1)  
%         lflag(it)=-1;
%     else
%         lflag(it)=0;
%     end
%     
%     if test(it,end)>test(it-1,end)
%         rflag(it)=1;
%     elseif test(it,end)<test(it-1,end) 
%         rflag(it)=-1;
%     else
%         rflag(it)=0;
%     end
% end
% nlp=0;
% nrp=0;
% for it=2:length(t)
%     if lflag(it)==-1 && lflag(it-1)==1
%         nlp=nlp+1;
%         lpeak(nlp)=t(it);
%     end
%     if rflag(it)==-1 && rflag(it-1)==1
%         nrp=nrp+1;
%         rpeak(nrp)=t(it);
%     end
% end

%right
%original detection of peaks
[pksr,locr,wr,pr]=findpeaks(test(:,end));%pks=peak value, loc=temporal index of peak, w=width of peak,p= prominence of peak
                                                            %use t(loc) for actual peak times
minpr=0;
sortpr=sort(pr,'descend');
for i=1:length(pr)-1
    if sortpr(i)/5>sortpr(i+1)
        minpr=sortpr(i);
        break;
    end
end
%new detection
[npksr,nlocr,nwr,npr]=findpeaks(test(:,end),'MinPeakProminence',minpr);%pks=peak value, loc=temporal index of peak, w=width of peak,p= prominence of peak
                                                            %use t(loc) for actual peak times
ntlocr=t(nlocr);

%left
[pksl,locl,wl,pl]=findpeaks(test(:,1));%pks=peak value, loc=temporal index of peak, w=width of peak,p= prominence of peak
                                                            %use t(loc) for actual peak times
minpl=0;
sortpl=sort(pl,'descend');
for i=1:length(pl)-1
    if sortpl(i)/5>sortpl(i+1)
        minpl=sortpl(i);
        break;
    end
end
%new detection
[npksl,nlocl,nwl,npl]=findpeaks(test(:,1),'MinPeakProminence',minpl);%pks=peak value, loc=temporal index of peak, w=width of peak,p= prominence of peak
                                                            %use t(loc) for actual peak times
ntlocl=t(nlocl);
                                                         

