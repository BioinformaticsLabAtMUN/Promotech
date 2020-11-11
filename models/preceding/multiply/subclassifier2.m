%the second one: identify sigma24 promoters
function [Name3,testD3,sample_label3,label3] = subclassifier2(Name2,testD2,sample_label2,label2)
mm=size(sample_label2,1);
sample_label3=sample_label2;
if all(label2==1)
    testD3=1;
    Name3=1;
    label3=1;
    for i=1:mm
        sample_label3(i,3)=0;
    end
    return
end
load('allfeatureL2_2.mat');
load F-scorezhiL2_2.mat
TrainFeatureVector=allfeatureL2_2;
TrainLabel=yapp;
PPT21 = KNN2(Name2,testD2);
PPT22 = BPB2(Name2,testD2);
PPT23 = DNC2(Name2,testD2);
PPT25 = DAC2(Name2,testD2);
n2=length(testD2);
xtest=[PPT22,PPT21,PPT25,PPT23];%特征加一块
ytest=zeros(n2,1);
% %---------------------------
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 2.8284 -g 2 -w1 1 -w-1 1');
[predict_label2, accuracy, dec_values]=svmpredict(ytest,xtest,model);

k=0;
for m=1:mm
    if (sample_label2(m,2)==-1)
        k=k+1;
        sample_label3(m,3) = predict_label2(k,1);
    else
        sample_label3(m,3)=0;
    end
end
%若都预测为sigma24就退出程序，输出结果
label3=zeros(n2,1);
j=0;
for m=1:n2
    if predict_label2(m)==-1
        j=j+1;
        testD3{1,j}=testD2{1,m};
        Name3{1,j}=Name2{1,m};
        label3(m) = -1;
    else
        label3(m) = 1;
    end
end 

if all(label3==1)
    testD3=1;
    Name3=1;
end
     
%     fidout=fopen('results.txt','w+');
%     for i=1:Total_number_fragments
%         a=find(sample_label3(i,:)==1);
%         if a==1
%             fprintf(fidout,'The query sequence (%s) is a non-promoter sequence\r\n',Name{1,i});
%         elseif a==2
%             fprintf(fidout,'The query sequence (%s) is a sigma70-promoter sequence\r\n',Name{1,i});
%         else
%             fprintf(fidout,'The query sequence (%s) is a sigma24-promoter sequence\r\n',Name{1,i});
%         end
%     end
%     fclose(fidout);
%     return
% end

%save subclassifier2 Name Total_number_fragments Name3 testD3
%the second one: identify sigma32 promoters
%[sample_label4] = subclassifier3(output,sample_label3);



