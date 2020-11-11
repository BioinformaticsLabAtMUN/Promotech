clc

root_path = "../../data/bacteria_1_1/";
bacteria  = ["MYCOBACTER" "CLOSTRIDIUM" "RHODOBACTER_1" "RHODOBACTER_2" "BACILLUS" ];

import java.util.ArrayList;
TPs  =  ArrayList();
FNs  =  ArrayList();
FPs  =  ArrayList();
TNs  =  ArrayList();
PREs =  ArrayList();
RECs =  ArrayList();
TPRs =  ArrayList();
FPRs =  ArrayList();
ACCs =  ArrayList();
F1s  =  ArrayList();
MCCs =  ArrayList();

PROBS = containers.Map;

disp(["RUNNING ", bacteria ]);

for b_i = 1: length(bacteria)
    tic
    current_bacteria        = bacteria(b_i);
    bacteria_output_dir     = strcat('data/', current_bacteria);
    bacteria_output_file_X  = strcat(bacteria_output_dir, "/X.mat");
    bacteria_output_file_Y  = strcat(bacteria_output_dir, "/Y.mat");
    bacteria_output_results      = strcat(bacteria_output_dir, "/RESULTS.csv");
    bacteria_output_file_SEQS    = strcat(bacteria_output_dir, "/SEQS.mat");
    bacteria_output_file_LABELS  = strcat(bacteria_output_dir, "/LABELS.mat");
    bacteria_output_file_CHROMS  = strcat(bacteria_output_dir, "/CHROMS.mat");
    load(  bacteria_output_file_X  );
    load(  bacteria_output_file_Y  );
    load(  bacteria_output_file_SEQS    ); 
    load(  bacteria_output_file_LABELS  );
    load(  bacteria_output_file_CHROMS  );  
    fprintf( "___________________LOAD DATA FOR %s___________________\n\nX: %s\nY: %s",current_bacteria, bacteria_output_file_X, bacteria_output_file_Y );
    fprintf('\nSHAPES: \n\nCHROM SHAPE: (%s) \nSEQ SHAPE: (%s) \nLABELS: (%s) \n', num2str(size(CHROMS)), num2str(size(SEQS)), num2str(size(LABELS)));
    fprintf('\nSAMPLES: \n\nCHROM: %s \nSEQ: %s \nLABEL: %d \n', CHROMS(1), SEQS(1), LABELS(1) );        
    fprintf('\nX: (%s) \nY: (%s) \n', num2str(size(X)), num2str(size(Y)));
    
    fprintf( "___________________LOAD MODEL___________________\n\n" );
    load('model.mat')
    disp(model)
    
    fprintf( "___________________PREDICT___________________\n\n" );
    [Y_PRED, ACC, PROB]=svmpredict(Y,X,model, '-b 1');
    Y_PRED(Y_PRED==-1)=0;
    
    fprintf("CHROMS SHAPE: (%s) SAMPLE: %s          \n"     , num2str(size(CHROMS)), num2str(CHROMS(1)) );
    fprintf("SEQS   SHAPE: (%s) SAMPLE: %s          \n"     , num2str(size(SEQS))  , num2str(SEQS(1))   );
    fprintf("X      SHAPE: (%s) SAMPLE: %s          \n"     , num2str(size(X))     , num2str(X(1)) );
    fprintf("Y      SHAPE: (%s) SAMPLE: %s          \n"     , num2str(size(Y))     , num2str(Y(1)) );
    fprintf("Y_PRED     SHAPE: (%s) SAMPLE: %s      \n"     , num2str(size(Y_PRED)) , num2str(Y_PRED(1:10)) );
    fprintf("ACC        SHAPE: (%s) SAMPLE: %s      \n"     , num2str(size(ACC))      , num2str(ACC)      );
    fprintf("DEC VALUES SHAPE: (%s) dec_values:   \n"     , num2str(size(PROB))     );
    disp(PROB(1:10) );
    
    PROBS(current_bacteria) = dec_values;
    fprintf("\nMIN %s MAX, %s\n\n", num2str(min(PROB)) , num2str(max(PROB)) );
    
    fprintf( "___________________SAVING RESULTS TO %s___________________\n\n", bacteria_output_results );
    T_RESULTS = table( CHROMS.', SEQS.', Y, Y_PRED, PROB, 'VariableNames',{ 'CHROM','SEQ', 'Y', 'Y_PRED', 'PROB' } );
    disp(head(T_RESULTS));
    writetable( T_RESULTS, bacteria_output_results ,'Delimiter',',');
    
    C = confusionmat(Y,Y_PRED, 'Order',[0 1]);
    disp("CONFUSION MATRIX. Y AXIS = TRUE, X AXIS = PREDICTED");
    disp(C)
    
    TP =  C(1,1);	
    FN =  C(1,2); 
	FP =  C(2,1); 	
    TN =  C(2,2); 
    
    fprintf("\n\nTP: %s \nFN: %s \nFP: %s \nTN: %s \n\n", num2str(TP), num2str(FN), num2str(FP), num2str(TN) );
    
    PRE = TP / (TP + FP); %precision or positive predictive value (PPV)
    REC = TP / (TP + FN); %sensitivity, recall, hit rate, or true positive rate (TPR)
    TPR = REC;
    FPR = FP / (FP + TN);
    ACC = ( FP + TN ) / ( TP + TN + FP + FN ) ;
    F1  = (2 * TP) / ( (2*TP) + FP + FN );
    MCC = ((TP * TN) - (FP * FN))/sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) );
    
    
    TPs.add( TP );
    FNs.add( FN );
    FPs.add( FP );
    TNs.add( TN );
    PREs.add( PRE );
    RECs.add( REC );
    TPRs.add( TPR );
    FPRs.add( FPR );
    ACCs.add( ACC );
    F1s.add( F1 );
    MCCs.add( MCC );

    fprintf("\n\nPRE: %s \nREC: %s \nTPR: %s \nFPR: %s \nACC: %s \nF1 :  %s \nMCC: %s \n\n", num2str(PRE), num2str(REC), num2str(TPR), num2str(FPR), num2str(ACC) , num2str(F1), num2str(MCC) );

    
    %confusionchart(C);
    %[ap, recall, precision] = evaluateDetectionPrecision(Y_PRED, Y);
    %fprintf("\nPRE: %s \nREC: %s \nAP: %s \n\n", num2str(precision), num2str(recall), num2str(ap) );
    %fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);
end

TPs  =  cell2mat(cell(toArray(TPs)));
FNs  =  cell2mat(cell(toArray(FNs)));
FPs  =  cell2mat(cell(toArray(FPs)));
TNs  =  cell2mat(cell(toArray(TNs)));
PREs =  cell2mat(cell(toArray(PREs)));
RECs =  cell2mat(cell(toArray(RECs)));
TPRs =  cell2mat(cell(toArray(TPRs)));
FPRs =  cell2mat(cell(toArray(FPRs)));
ACCs =  cell2mat(cell(toArray(ACCs)));
F1s  =  cell2mat(cell(toArray(F1s)));
MCCs =  cell2mat(cell(toArray(MCCs)));

fprintf( "___________________SUMMARY___________________\n\n" );
T_SUMMARY = table(bacteria.', TPs, FNs, FPs, TNs, PREs, RECs, TPRs, FPRs, ACCs, F1s, MCCs, 'VariableNames',{ 'BACTERIA', 'TP', 'FN', 'FP', 'TN', 'PRE', 'REC', 'TPR', 'FPR', 'ACC', 'F1', 'MCC' } );
disp(T_SUMMARY);
writetable(T_SUMMARY,'data/RESULTS.csv','Delimiter',',')

fprintf( "___________________PROBS___________________\n\n" );
disp(keys(PROBS));
disp(values(PROBS));




% sample_label1=zeros(Total_number_fragments,6);
% j=0;
% for m=1:Total_number_fragments   
%     if predict_label(m)==1
%         j=j+1;
%         testD1{1,j}=testD{1,m};
%         Name1{1,j}=Name{1,m};
%         sample_label1(m,1)=0;
%         label1(m)=0;
%     else
%         sample_label1(m,1)=1;
%         label1(m)=1;
%     end
% end
% if all(label1==1)
%     testD1=1;
%     Name1=1;
% end
% %further verify which of the six types the identified promoter belongs to
% %the first one: identify sigma70 promoters
% [Name2,testD2,sample_label2,label2] = subclassifier1(Name1,testD1,sample_label1,label1);
% %the second one: identify sigma24 promoters
% [Name3,testD3,sample_label3,label3] = subclassifier2(Name2,testD2,sample_label2,label2);
% %the third one: identify sigma32 promoters
% [Name4,testD4,sample_label4,label4] = subclassifier3(Name3,testD3,sample_label3,label3);
% %the fourth one: identify sigma38 promoters
% [Name5,testD5,sample_label5,label5] = subclassifier4(Name4,testD4,sample_label4,label4);
% %the fifth one: identify sigma28 promoters
% [sample_label6] = subclassifier5(Name5,testD5,sample_label5,label5);
% fidout=fopen('results.txt','w');    %Please change here to specify the output results file.
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




