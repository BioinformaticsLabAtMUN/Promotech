function PPT3 = DNC(Name,testD)
%˫�˸��Ậ��
load F-scorezhiL1.mat F3
hn=Name; 
negative=testD;
Nn=length(hn);%number of negative samples
for i=1:Nn
    Str=negative{1,i};
    Str=char(Str);
    Str=upper(Str);%�����������е�Сд����һ�ɻ��ɴ�д
    negative{1,i}=Str;
end
AA='ACGT';
L=length(testD{1,1});%ÿ�����еĳ���
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
for i=1:9
    k=F3(2,i);
    PPT3(:,i)=testDNC(:,k);%����������������������F-score������������
end
