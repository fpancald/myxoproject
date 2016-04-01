cd teststat2d/diffusiontest
load('DD1.mat');
avgDD=DD;
minDD=DD;
maxDD=DD;
dim1=size(DD,1);
dim2=size(DD,2);
for i=2:100
    filename=strcat(strcat('DD',num2str(i)),'.mat');
    load(filename);
    avgDD=avgDD+DD;
    for i1=1:dim1
        for i2=1:dim2
            if minDD(i1,i2)>DD(i1,i2)
                minDD(i1,i2)=DD(i1,i2);
            end;
            if maxDD(i1,i2)<DD(i1,i2)
                maxDD(i1,i2)=DD(i1,i2);
            end;
        end
    end
    clear DD;
end
avgDD=avgDD/100;
cd ../..