function PPT2 = BPB(Name,testD)
load F-scorezhiL1.mat F2
[hp,positive]=fastaread('positive2860.txt');%trainpositive.txt
Np=length(hp);
[hn,negative]=fastaread('negative2860.txt');
Nn=length(hn);
x1=Name;
WNp=length(x1);%number of positive samples
%每个蛋白质肽段的长度
for i=1:WNp
    Str=testD{1,i};
    Str=char(Str);
    Str=upper(Str);%将测试样本中的小写序列一律换成大写
    testD{1,i}=Str;
    %testD1(i,:)=Str;
end
for i=1:Np
    Str=positive{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive{1,i}=Str;
end
for i=1:Nn
    Str=negative{1,i};
    Str=char(Str);
    negative{1,i}=Str;
end
%—————————————————————————
%——————————————————————————
AA='ACGT';
%记录每个氨基酸在每个位置出现的次数
M=length(negative{1,1});
F11=zeros(4,M);%记录每个氨基酸在每个位置出现的次数
F12=zeros(4,M);
for m=1:Np
    for j=1:M
        t=positive{1,m}(j);
        k=strfind(AA,t);
        F11(k,j)=F11(k,j)+1;
    end
end
%
for m=1:Nn
    for j=1:M
        t=negative{1,m}(j);
        k=strfind(AA,t);
        F12(k,j)=F12(k,j)+1;
    end
end
F11= F11/Np;
F12= F12/Nn;

testBPB=zeros(WNp,2*M);
for m=1:WNp
    l=M;
    for j=1:l
        %%%%%%%%%%%%%%%%%%%%%%%%
        if(j < length(testD{1,m})  )
%             fprintf("%s %s %s \n", num2str(m),num2str(j), num2str(length(testD{1,m}) )  )  
%             fprintf("%s\n", num2str(testD{1,m}(j)) )
            t=testD{1,m}(j);
            k=strfind(AA,t);
            testBPB(m,j)=F11(k,j);
            testBPB(m,l+j)=F12(k,j);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
%         t=testD{1,m}(j);
%         k=strfind(AA,t);
%         testBPB(m,j)=F11(k,j);
%         testBPB(m,l+j)=F12(k,j);
     end
end
for i=1:130
    k=F2(2,i);
    PPT2(:,i)=testBPB(:,k);%将正负样本的所有特征按F-score分数进行排序
end







