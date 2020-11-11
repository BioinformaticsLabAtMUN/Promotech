%the first one: identify sigma70 promoters
function [Name2, testD2,sample_label2,label2] = subclassifier1(Name1, testD1, sample_label1,label1)
mm=size(sample_label1,1);
sample_label2=sample_label1;
if all(label1==1)
    testD2=1;
    Name2=1;
    label2=1;
    for i=1:mm
        sample_label2(i,2)=0;
    end
    return
end
load('allfeatureL2_1.mat');
load F-scorezhiL2_1.mat
TrainFeatureVector=allfeatureL2_1;
TrainLabel=yapp;
PPT11 = KNN1(Name1,testD1);
PPT12 = BPB1(Name1,testD1);
PPT15 = DAC1(Name1,testD1);
n1=length(testD1);
xtest=[PPT11,PPT12,PPT15];%特征加一块
ytest=zeros(n1,1);
% %---------------------------
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 1.4142 -g 0.01121 -w1 1 -w-1 1');
[predict_label1, accuracy, dec_values]=svmpredict(ytest,xtest,model);

k=0;
for m=1:mm
    if (sample_label1(m,1)==0)
        k=k+1;
        sample_label2(m,2) = predict_label1(k,1);
    else
        sample_label2(m,2)=0;
    end
end
%若都预测为sigma70就退出程序，输出结果

label2=zeros(n1,1);
j=0;
for m=1:n1
    if predict_label1(m)==-1
        j=j+1;
        testD2{1,j}=testD1{1,m};
        Name2{1,j}=Name1{1,m};
        label2(m)=-1;
    else
        label2(m)=1;
    end
end 

if all(label2==1)
    testD2=1;
    Name2=1;
end
%     fidout=fopen('results.txt','w+');
%     for i=1:Total_number_fragments
%         a=find(sample_label2(i,:)==1);
%         if a==1
%             fprintf(fidout,'The query sequence (%s) is a non-promoter sequence\r\n',Name{1,i});
%         else
%             fprintf(fidout,'The query sequence (%s) is a sigma70-promoter sequence\r\n',Name{1,i});
%         end
%     end
%     fclose(fidout);
%     return
% end

%save subclassifier1 Name Total_number_fragments Name2 testD2
%the second one: identify sigma24 promoters
%[sample_label3] = subclassifier2(sample_label2);
        



