%the fifth one: identify sigma28 promoters
function [sample_label6] = subclassifier5(Name5,testD5,sample_label5,label5)
mm=size(sample_label5,1);
sample_label6=sample_label5;
if all(label5==1)
    for i=1:mm
        sample_label6(i,6)=0;
    end
    return
end
load('allfeatureL2_5.mat');
load F-scorezhiL2_5.mat
TrainFeatureVector=allfeatureL2_5;
TrainLabel=yapp;
PPT51 = KNN5(Name5,testD5);
PPT52 = BPB5(Name5,testD5);
PPT53 = DNC5(Name5,testD5);
PPT55 = DAC5(Name5,testD5);
n5=length(testD5);
xtest=[PPT52,PPT51,PPT53,PPT55];%ÌØÕ÷¼ÓÒ»¿é
ytest=zeros(n5,1);
% %---------------------------
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 1.4142 -g 2 -w1 1 -w-1 1');
[predict_label5, accuracy, dec_values]=svmpredict(ytest,xtest,model);

k=0;
for m=1:mm
    if (sample_label5(m,5)==-1)
        k=k+1;
        sample_label6(m,6) = predict_label5(k,1);
    else
        sample_label6(m,6)=0;
    end
end
% fidout=fopen('results.txt','w+');
% for i=1:Total_number_fragments
%     a=find(sample_label6(i,:)==1);
%     if a==1
%        fprintf(fidout,'The query sequence (%s) is a non-promoter sequence\r\n',Name{1,i});
%     elseif a==2
%        fprintf(fidout,'The query sequence (%s) is a sigma70-promoter sequence\r\n',Name{1,i});
%     elseif a==3
%        fprintf(fidout,'The query sequence (%s) is a sigma24-promoter sequence\r\n',Name{1,i});
%     elseif a==4
%        fprintf(fidout,'The query sequence (%s) is a sigma32-promoter sequence\r\n',Name{1,i});
%     elseif a==5
%        fprintf(fidout,'The query sequence (%s) is a sigma38-promoter sequence\r\n',Name{1,i});
%     elseif a==6
%        fprintf(fidout,'The query sequence (%s) is a sigma28-promoter sequence\r\n',Name{1,i});
%     else
%        fprintf(fidout,'The query sequence (%s) is a sigma54-promoter sequence\n',Name{1,i});
%     end           
% end
% fclose(fidout);



