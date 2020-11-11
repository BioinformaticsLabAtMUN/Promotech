function PPT1 =KNN(Name,testD)
%KNN 
load F-scorezhiL1.mat F1
[hp,positive]=fastaread('positive2860.txt');
[hn,negative]=fastaread('negative2860.txt');
Np=length(hp);%number of positive samples
Nn=length(hn);%number of negative samples
x1=Name;
WNp=length(x1);%number of positive samples
%The length of each protein peptide

disp("Change all lowercase sequences in the test sample to uppercase")
for i=1:WNp
    Str=testD{1,i};
    Str=char(Str);
    Str=upper(Str);
    testD{1,i}=Str;
    testD1(i,:)=Str;
end
disp("Replace all lowercase sequences in positive samples with uppercase")
for i=1:Np
    Str=positive{1,i};
    Str=char(Str);
    Str=upper(Str); 
    positive{1,i}=Str;
    ppp(i,:)=Str;
end
disp("Replace all lowercase sequences in negative samples with uppercase")
for i=1:Nn
    Str=negative{1,i};
    Str=char(Str);
    negative{1,i}=Str;
    nnn(i,:)=Str;
end
M1=length(testD{1,1});
xp=[ppp;nnn];
M2=length(xp(:,1));
ss1=[1,M2];
o=0;
disp("Starting KNN")
f = waitbar(0, sprintf( "STARTING KNN 0/200") );
for k=10:10:200
    o=o+1;
    for j=1:WNp
        for i=1:M2
            a=0;
            for ii=1:M1
                if strcmp(testD1(j,ii),xp(i,ii))
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
        waitbar(k/200,f,sprintf( "%s/%s - %s/%s", num2str(k), num2str(200), num2str(j), num2str(WNp) ))
    end 
end
close(f)
for i=1:15
    k=F1(2,i);
    PPT1(:,i)=w1(:,k); %Sort all features of positive and negative samples by F-score score
end


