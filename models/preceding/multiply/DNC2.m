function PPT23 = DNC2(Name2,testD2)
%双核甘酸含量
load F-scorezhiL2_2.mat F3
hn=Name2; 
negative=testD2;
Nn=length(hn);%number of negative samples
for i=1:Nn
    Str=negative{1,i};
    Str=char(Str);
    Str=upper(Str);%将测试样本中的小写序列一律换成大写
    negative{1,i}=Str;
end
AA='ACGT';
L=length(testD2{1,1});%每个序列的长度
testDNC=zeros(Nn,16);
for i=1:Nn
    Di_nucleotides=zeros(4,4);
    Peptide=negative{1,i};
    for j=1:L-1
        s=Peptide(j);
        r=strfind(AA,s);
        t=Peptide(j+1);
        l=strfind(AA,t);
        Di_nucleotides(r,l)=Di_nucleotides(r,l)+1;
    end    
    Di_nucleotides=Di_nucleotides/(L-1);
    A=Di_nucleotides';
    B=A(:);
    testDNC(i,:)=B';
end
for i=1:10
    k=F3(2,i);
    PPT23(:,i)=testDNC(:,k);%将正负样本的所有特征按F-score分数进行排序
end

