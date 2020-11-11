function PPT51 =KNN5(Name5,testD5)
load F-scorezhiL2_5.mat F1
[hp2,positive2]=fastaread('sigma28promoter.txt');
Np2=length(hp2);
[hp5,positive5]=fastaread('sigma54promoter.txt');
Np5=length(hp5);

for i=1:Np2
    Str=positive2{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive2{1,i}=Str;
    ppp2(i,:)=Str;
end
for i=1:Np5
    Str=positive5{1,i};
    Str=char(Str);
    Str=upper(Str);%将正样本中的小写序列一律换成大写
    positive5{1,i}=Str;
    ppp5(i,:)=Str;
end

ppp=ppp2;
nnn=ppp5;
Np=length(ppp);
Nn=length(nnn);

x1=Name5;
WNp=length(x1);%number of positive samples

for i=1:WNp
    Str=testD5{1,i};
    Str=char(Str);
    Str=upper(Str);%将测试样本中的小写序列一律换成大写
    testD5{1,i}=Str;
    testD55(i,:)=Str;
end
M1=length(testD5{1,1});
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
                if strcmp(testD55(j,ii),xp(i,ii))
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
for i=1:3
    k=F1(2,i);
    PPT51(:,i)=w1(:,k);%将正负样本的所有特征按F-score分数进行排序
end
