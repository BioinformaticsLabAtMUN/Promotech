clc
%%%%function  PromoterPredweb
%%%%promoter_web-predict 
%%%%%For the first layer, we selected the feature combination KNN(15)+BPB(130)+DNC(9)+MNC(1)+DAC(10)

load('allfeatureL1.mat')
load('F-scorezhiL1.mat')

TrainFeatureVector=allfeatureL1;
TrainLabel=yapp;

[head,seq]=fastaread('sample.fasta');      %Please change here to specify the input file to PromoterPred2L
seq=cellstr(seq);
head=cellstr(head);
SeqNumber=length(seq); %The number of sequences
%--------------------------------------------
%%identify the possible promoter sequence
%% Ectract DNA fragments with a length of 81
t=0;

disp(["CUTTING INTO : ", SeqNumber])

f = waitbar(0,'Please wait...');
for i=1:SeqNumber   
    Se=seq{1,i};%The i-th sample sequence
    L=length(Se);%The length of this sequence
    if L>=81
        Number_Length_81=L-81+1;%The number of fragments of length 81 in the i-th sample
        for j=1:Number_Length_81
            fragment_j=Se(j:j+80);%��ȡ��81Ƭ�Σ�
            Name{1,t+j}=char(head(i));
            testD{1,t+j}=(fragment_j);
        end   
        t=t+Number_Length_81;
    end
    waitbar(i/SeqNumber,f,sprintf(Se))
end
close(f)

Total_number_fragments=length(Name);%һ����ȡ������

disp(["# OF FRAGMENTS : ", Total_number_fragments])
AA='ACGT';
V=81;%the length of each sample
%%%%%%%%%%%%%%%%%%%%%%%
f = waitbar(0,'RUNNING KNN');
PPT1 = KNN(Name,testD);
waitbar(0.15,f,sprintf("RUNNING BPB"))
PPT2 = BPB(Name,testD);
waitbar(0.25,f,sprintf("RUNNING DNC"))
PPT3 = DNC(Name,testD);
waitbar(0.5,f,sprintf("RUNNING MNC"))
PPT4 = MNC(Name,testD);
waitbar(0.75,f,sprintf("RUNNING DAC"))
PPT5 = DAC(Name,testD);
waitbar(1,f,sprintf("FINISHED"))
close(f)


% 60 FRAGMENTS
xtest=[PPT1,PPT2,PPT3,PPT4,PPT5]; %������һ��
ytest=zeros(Total_number_fragments,1);  

load('model.mat')
disp("PREDICT")
%model=fitcsvm(TrainLabel, TrainFeatureVector );
[predict_label, accuracy, dec_values]=svmpredict(ytest,xtest,model);

sample_label1=zeros(Total_number_fragments,6);

j=0;
for m=1:Total_number_fragments   
    if predict_label(m)==1
        j=j+1;
        testD1{1,j}=testD{1,m};
        Name1{1,j}=Name{1,m};
        sample_label1(m,1)=0;
        label1(m)=0;
    else
        sample_label1(m,1)=1;
        label1(m)=1;
    end
end

if all(label1==1)
    testD1=1;
    Name1=1;
end

%%
%further verify which of the six types the identified promoter belongs to
%the first one: identify sigma70 promoters
[Name2,testD2,sample_label2,label2] = subclassifier1(Name1,testD1,sample_label1,label1);
%the second one: identify sigma24 promoters
[Name3,testD3,sample_label3,label3] = subclassifier2(Name2,testD2,sample_label2,label2);
%the third one: identify sigma32 promoters
[Name4,testD4,sample_label4,label4] = subclassifier3(Name3,testD3,sample_label3,label3);
%the fourth one: identify sigma38 promoters
[Name5,testD5,sample_label5,label5] = subclassifier4(Name4,testD4,sample_label4,label4);
%the fifth one: identify sigma28 promoters
[sample_label6] = subclassifier5(Name5,testD5,sample_label5,label5);

fidout=fopen('results.txt','w');    %Please change here to specify the output results file.

for i=1:Total_number_fragments
    a=find(sample_label6(i,:)==1);
    if a==1
       fprintf(fidout,'The query sequence (%s) is a non-promoter sequence\r\n',Name{1,i});
    elseif a==2
       fprintf(fidout,'The query sequence (%s) is a sigma70-promoter sequence\r\n',Name{1,i});
    elseif a==3
       fprintf(fidout,'The query sequence (%s) is a sigma24-promoter sequence\r\n',Name{1,i});
    elseif a==4
       fprintf(fidout,'The query sequence (%s) is a sigma32-promoter sequence\r\n',Name{1,i});
    elseif a==5
       fprintf(fidout,'The query sequence (%s) is a sigma38-promoter sequence\r\n',Name{1,i});
    elseif a==6
       fprintf(fidout,'The query sequence (%s) is a sigma28-promoter sequence\r\n',Name{1,i});
    else
       fprintf(fidout,'The query sequence (%s) is a sigma54-promoter sequence\n',Name{1,i});
    end           
end
fclose(fidout);










