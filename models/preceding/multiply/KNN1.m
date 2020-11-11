function PPT11 =KNN1(Name1,testD1)
%KNN 
load F-scorezhiL2_1.mat F1
[hp1,positive1]=fastaread('sigma24promoter.txt');%trainpositive.txt
Np1=length(hp1);
[hp2,positive2]=fastaread('sigma28promoter.txt');
Np2=length(hp2);
[hp3,positive3]=fastaread('sigma32promoter.txt');
Np3=length(hp3);
[hp4,positive4]=fastaread('sigma38promoter.txt');
Np4=length(hp4);
[hp5,positive5]=fastaread('sigma54promoter.txt');
Np5=length(hp5);
[hp6,positive6]=fastaread('sigma70promoter.txt');
Np6=length(hp6);

for i=1:Np1
    Str=positive1{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive1{1,i}=Str;
    ppp1(i,:)=Str;
end
for i=1:Np2
    Str=positive2{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive2{1,i}=Str;
    ppp2(i,:)=Str;
end
for i=1:Np3
    Str=positive3{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive3{1,i}=Str;
    ppp3(i,:)=Str;
end
for i=1:Np4
    Str=positive4{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive4{1,i}=Str;
    ppp4(i,:)=Str;
end
for i=1:Np5
    Str=positive5{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive5{1,i}=Str;
    ppp5(i,:)=Str;
end
for i=1:Np6
    Str=positive6{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive6{1,i}=Str;
    ppp6(i,:)=Str;
end
ppp=ppp6;
nnn=[ppp1;ppp2;ppp3;ppp4;ppp5];
Np=length(ppp);
Nn=length(nnn);

x1=Name1;
WNp=length(x1);%number of positive samples
%每个蛋白质肽段的长度
for i=1:WNp
    Str=testD1{1,i};
    Str=char(Str);
    Str=upper(Str);%将测试样本中的小写序列一律换成大写
    testD1{1,i}=Str;
    testD11(i,:)=Str;
end
p1=length(positive1{1,1});
M1=length(testD1{1,1});
xp=[ppp;nnn];
M2=length(xp(:,1));
ss1=[1,M2];
o=0;
for k=10:10:200
    o=o+1;
    for j=1:WNp
        for i=1:M2
            a=0;
            for ii=1:M1
                if strcmp(testD11(j,ii),xp(i,ii))
                    a=a+2;
                else
                    a=a-1;
                end
            end
            ss1(1,i)=a;
        end
        s1=-ss1;
        [q1,order1]=sort(s1);
        h1=0;
        for i=1:k+1
            if order1(i)<=Np
                h1=h1+1;
            end
        end
        w1(j,o)=h1/(k+1);
    end 
end
for i=1:15
    k=F1(2,i);
    PPT11(:,i)=w1(:,k);%将正负样本的所有特征按F-score分数进行排序
end
