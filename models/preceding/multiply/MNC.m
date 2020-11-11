function PPT4 = MNC(Name,testD)
%单核苷酸的频率
load F-scorezhiL1.mat F4
hp=Name;
positive=testD;
Np=length(hp);%number of positive samples

for i=1:Np
    Str=positive{1,i};
    Str=char(Str);
    Str=upper(Str);%将测试样本中的小写序列一律换成大写
    positive{1,i}=Str;
end
AA='ACGT';
L=length(testD{1,1});%每个肽段的长度
testMNC=zeros(Np,4);
for i=1:Np
    Peptide=positive{1,i};
    for j=1:L
        s=Peptide(j);
        k=strfind(AA,s);
        testMNC(i,k)= testMNC(i,k)+1;
    end
end
testMNC=testMNC/L;
for i=1:1
    k=F4(2,i);
    PPT4(:,i)=testMNC(:,k);%将正负样本的所有特征按F-score分数进行排序
end
