s=size(dc);
DC=[];
for it=1:s(1)
    DCit=[];
    for i=1:s(2)
        DCit=[DCit;dc{it,i}];
    end
    [Y,I]=sort(DCit(:,1));
    DCit=DCit(I,:); %use the column indices from sort() to sort all columns of A.
    DC{it}=DCit;
end
% DC

% DC

X=[];
D=[];
for it=1:s(1)
    M=DC{it};
    x=[];
    d=[];
    for i=1:s(2)
        if any(abs(x-M(i,1))<=1e-1)
            d(end)=d(end)+M(i,2);
        else
            x=[x M(i,1)];
            d=[d M(i,2)];
        end
    end
    X{it}=x;
    D{it}=d;
end
