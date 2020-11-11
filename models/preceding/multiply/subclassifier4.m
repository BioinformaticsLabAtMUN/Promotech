%the fourth one: identify sigma38 promoters
function [Name5,testD5,sample_label5,label5] = subclassifier4(Name4,testD4,sample_label4,label4)
mm=size(sample_label4,1);
sample_label5=sample_label4;
if all(label4==1)
    testD5=1;
    Name5=1;
    label5=1;
    for i=1:mm
        sample_label5(i,5)=0;
    end
    return
end
load('allfeatureL2_4.mat');
load F-scorezhiL2_4.mat
TrainFeatureVector=allfeatureL2_4;
TrainLabel=yapp;
PPT41 = KNN4(Name4,testD4);
PPT42 = BPB4(Name4,testD4);
n4=length(testD4);
xtest=[PPT41,PPT42];%特征加一块
ytest=zeros(n4,1);
% %---------------------------
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 32 -g 0.25 -w1 1 -w-1 1');
[predict_label4, accuracy, dec_values]=svmpredict(ytest,xtest,model);

k=0;
for m=1:mm
    if (sample_label4(m,4)==-1)
        k=k+1;
        sample_label5(m,5) = predict_label4(k,1);
    else
        sample_label5(m,5)=0;
    end
end
%若都预测为sigma38就退出程序，输出结果
label5=zeros(n4,1);
j=0;
for m=1:n4
    if predict_label4(m)==-1
        j=j+1;
        testD5{1,j}=testD4{1,m};
        Name5{1,j}=Name4{1,m};
        label5(m)=-1;
    else
        label5(m)=1;
    end 
end

if all(label5==1)
    testD5=1;
    Name5=1;
end
%     fidout=fopen('results.txt','w+');
%     for i=1:Total_number_fragments
%         a=find(sample_label5(i,:)==1);
%         if a==1
%             fprintf(fidout,'The query sequence (%s) is a non-promoter sequence\r\n',Name{1,i});
%         elseif a==2
%             fprintf(fidout,'The query sequence (%s) is a sigma70-promoter sequence\r\n',Name{1,i});
%         elseif a==3
%             fprintf(fidout,'The query sequence (%s) is a sigma24-promoter sequence\r\n',Name{1,i});
%         elseif a==4
%             fprintf(fidout,'The query sequence (%s) is a sigma32-promoter sequence\r\n',Name{1,i});
%         else
%             fprintf(fidout,'The query sequence (%s) is a sigma38-promoter sequence\r\n',Name{1,i});
%         end
%     end
%     fclose(fidout);
%     return
% end

%save subclassifier4 Name Total_number_fragments Name5 testD5
% %the second one: identify sigma28 promoters
% [sample_label6] = subclassifier5(output,sample_label5);




