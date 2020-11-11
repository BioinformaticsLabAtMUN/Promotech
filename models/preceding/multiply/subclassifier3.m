%the third one: identify sigma32 promoters
function [Name4,testD4,sample_label4,label4] = subclassifier3(Name3,testD3,sample_label3,label3)
mm=size(sample_label3,1);
sample_label4=sample_label3;
if all(label3==1)
    testD4=1;
    Name4=1;
    label4=1;
    for i=1:mm
        sample_label4(i,4)=0;
    end
    return
end
load('allfeatureL2_3.mat');
load F-scorezhiL2_3.mat
TrainFeatureVector=allfeatureL2_3;
TrainLabel=yapp;
PPT31 = KNN3(Name3,testD3);
PPT32 = BPB3(Name3,testD3);
PPT33 = DNC3(Name3,testD3);
n3=length(testD3);
xtest=[PPT32,PPT31,PPT33];%特征加一块
ytest=zeros(n3,1);
% %---------------------------
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 5.6569 -g 1 -w1 1 -w-1 1');
[predict_label3, accuracy, dec_values]=svmpredict(ytest,xtest,model);

k=0;
for m=1:mm
    if (sample_label3(m,3)==-1)
        k=k+1;
        sample_label4(m,4) = predict_label3(k,1);
    else
        sample_label4(m,4)=0;
    end
end
%若都预测为sigma32就退出程序，输出结果
label4=zeros(n3,1);
j=0;
for m=1:n3
    if predict_label3(m)==-1
        j=j+1;
        testD4{1,j}=testD3{1,m};
        Name4{1,j}=Name3{1,m};
        label4(m)=-1;
    else
        label4(m)=1;
    end    
end

 if all(label4==1)
     Name4=1;
     testD4=1;
 end
%     fidout=fopen('results.txt','w+');
%     for i=1:Total_number_fragments
%         a=find(sample_label4(i,:)==1);
%         if a==1
%             fprintf(fidout,'The query sequence (%s) is a non-promoter sequence\r\n',Name{1,i});
%         elseif a==2
%             fprintf(fidout,'The query sequence (%s) is a sigma70-promoter sequence\r\n',Name{1,i});
%         elseif a==3
%             fprintf(fidout,'The query sequence (%s) is a sigma24-promoter sequence\r\n',Name{1,i});
%         else
%             fprintf(fidout,'The query sequence (%s) is a sigma32-promoter sequence\r\n',Name{1,i});
%         end
%     end
%     fclose(fidout);
%     return
% end

% save subclassifier3 Name Total_number_fragments Name4 testD4
% %the second one: identify sigma38 promoters
% [sample_label5] = subclassifier4(output,sample_label4);


