function PPT55 = DAC5(Name5,testD5)
load X_norm.mat M
load F-scorezhiL2_2.mat F5
x1=Name5;
Np=length(x1);
for i=1:Np
    Str=testD5{1,i};
    Str=char(Str);
    Str=upper(Str);%将测试样本中的小写序列一律换成大写
    testD5{1,i}=Str;
end
L=length(testD5{1,1});%每个肽段的长度
AA='ACGT';
%%%%正样本
for i=1:Np
    Peptide=testD5{1,i};
    F1=zeros(6,1);
    for j=1:L-1
        t1=Peptide(j);
        k1=strfind(AA,t1);
        t2=Peptide(j+1);
        k2=strfind(AA,t2);
        F1=F1+M(:,4*(k1-1)+k2);
    end
    pp1(i,:)=F1'/(L-1);
end
lag=1;
for i=1:Np
    Value1=zeros(6,1);
    Peptide=testD5{1,i};
    Mean=pp1(i,:);
    for j=1:L-lag-1
        t1=Peptide(j);
        k1=strfind(AA,t1);
        t2=Peptide(j+1);
        k2=strfind(AA,t2);
        Index_k=4*(k1-1)+k2;
        t1_lag=Peptide(j+lag);
        k1_lag=strfind(AA,t1_lag);
        t2_lag=Peptide(j+lag+1);
        k2_lag=strfind(AA,t2_lag);
        Index_lag=4*(k1_lag-1)+k2_lag;
        Value1=Value1+(M(:,Index_k)-Mean').*(M(:,Index_lag)-Mean');  
    end
    m1(i,:)=Value1'/(L-lag-1);
end
lag=2;
for i=1:Np
    Value2=zeros(6,1);
    Peptide=testD5{1,i};
    Mean=pp1(i,:);
    for j=1:L-lag-1
        t1=Peptide(j);
        k1=strfind(AA,t1);
        t2=Peptide(j+1);
        k2=strfind(AA,t2);
        Index_k=4*(k1-1)+k2;
        t1_lag=Peptide(j+lag);
        k1_lag=strfind(AA,t1_lag);
        t2_lag=Peptide(j+lag+1);
        k2_lag=strfind(AA,t2_lag);
        Index_lag=4*(k1_lag-1)+k2_lag;
        Value2=Value2+(M(:,Index_k)-Mean').*(M(:,Index_lag)-Mean');   
    end
    m2(i,:)=Value2'/(L-lag-1);
end
testDAC=[m1,m2];
for i=1:3
    k=F5(2,i);
    PPT55(:,i)=testDAC(:,k);%将正负样本的所有特征按F-score分数进行排序
end




    

















