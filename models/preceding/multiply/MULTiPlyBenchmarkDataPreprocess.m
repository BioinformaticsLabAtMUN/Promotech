clc % CLEAR COMMAND LINE

%STEPS
% 1. Generate Sequence and Labels Array Using
%     generating_data   = true;
%     preprocess_data   = false;
%     generate_KNN_data = false;
% 2. Preprocess KNN PPT1
%     generating_data = false;
%     preprocess_data = true;
%     generate_KNN_data = false;
% 3. Generate PPT1,2,3,4 and joined them as X & Y
%     generating_data = false;
%     preprocess_data = true;
%     generate_KNN_data = true;
    

PREPROCESS();

function PREPROCESS(bacteria_to_run)
    generating_data = false;
    preprocess_data = true;
    generate_KNN_data = false;
    bacteria = ["MYCOBACTER" "CLOSTRIDIUM" "RHODOBACTER_1" "RHODOBACTER_2" "BACILLUS" ];
%     bacteria = [ "CLOSTRIDIUM" ];
    root_path = "../../data/bacteria_1_1/";
    
    if( exist('bacteria_to_run', 'var') )
        disp("RUNNING ALL BACTERIA")
%     else
%         clear bacteria
%         bacteria  = [ bacteria_to_run ];
%         disp(["RUNNING ONLY", bacteria, size(bacteria)])
    end
    
    disp(["RUNNING ", bacteria ]);
    
    for b_i = 1: length(bacteria)
        current_bacteria = bacteria(b_i);

        p_path = strcat(root_path, current_bacteria, "/positive.fasta");
        n_path = strcat(root_path, current_bacteria, "/negative.fasta");


        fprintf("*********************************%s**********************************\nP: %s \nN: %s \n", current_bacteria, p_path, n_path);
        bacteria_output_dir = strcat('data/', current_bacteria);
        [status, msg, msgID] = mkdir(bacteria_output_dir);
        disp(msg);


        bacteria_output_file_X = strcat(bacteria_output_dir, "/X.mat");
        bacteria_output_file_Y = strcat(bacteria_output_dir, "/Y.mat");
        bacteria_output_file_SUMMARY = strcat(bacteria_output_dir, "/SUMMARY.mat");
        bacteria_output_file_SEQS    = strcat(bacteria_output_dir, "/SEQS.mat");
        bacteria_output_file_LABELS  = strcat(bacteria_output_dir, "/LABELS.mat");
        bacteria_output_file_CHROMS  = strcat(bacteria_output_dir, "/CHROMS.mat");

        if(generating_data == true)
            disp("___________________LOADING POSITIVE SAMPLES___________________")
            tic
            [p_head,p_seq] = fastaread(p_path); %'test_sample.txt'
            p_seq     = cellstr(p_seq);
            p_head    = cellstr(p_head);
            p_seq_len = length(p_seq);
            p_labels  = ones(1, p_seq_len);
            PT        = table( p_head.', p_seq.', p_labels.', 'VariableNames',{'CHROM','SEQ', 'LABEL'} );
            disp( head(PT, 3) );
            fprintf('# POSITIVE SEQS: %s\n', num2str(p_seq_len));
            fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);

            disp("___________________LOADING NEGATIVE SAMPLES___________________");
            tic
            [n_head,n_seq] = fastaread(n_path); %'negative2860.txt'
            n_seq     = cellstr(n_seq);
            n_head    = cellstr(n_head);
            n_seq_len = length(n_seq);
            n_labels  = zeros(1, n_seq_len);
            NT        = table( n_head.', n_seq.', n_labels.', 'VariableNames',{'CHROM','SEQ', 'LABEL'}  );
            % NT        = NT(1:p_seq_len, :); %MAKE IT BALANCED
            % n_seq_len = p_seq_len;          %UPDATE TABLE NUMBER OF ROWS
            disp( head(NT, 3) );
            fprintf('# NEGATIVE SEQS: %s\n', num2str(n_seq_len));
            fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);


            disp("___________________COMBINE SAMPLES___________________");
            tic
            WT = [PT; NT];
            both_size_margin_sample = 3;
            disp( WT( p_seq_len-both_size_margin_sample:p_seq_len+both_size_margin_sample, : ) );
            fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);


            disp("___________________SHUFFLE SAMPLES___________________");
            tic
            SHUFFLE_T = WT(randperm(size(WT, 1)), :);
            CHROMS    = SHUFFLE_T{:,'CHROM'}.';
            SEQS      = SHUFFLE_T{:,'SEQ'}.';
            LABELS    = SHUFFLE_T{:,'LABEL'}.';
            shuff_len = length(SEQS);
            disp( head(SHUFFLE_T) );
            fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);


            disp("___________________SLIDING WINDOW OF 81BP___________________");
            tic
            fprintf('\nCHROM SHAPE: (%s) \nSEQ SHAPE: (%s) \nLABELS: (%s) \n', num2str(size(CHROMS)), num2str(size(SEQS)), num2str(size(LABELS)));
            t=0;
            import java.util.ArrayList;
            NEW_SEQS   = ArrayList();
            NEW_LABELS = ArrayList();
            NEW_CHROMS = ArrayList();
            f = waitbar(0,'STARTTING SLIDING WINDOW');
            for i=1:shuff_len   
                Se = char(SEQS(1,i));                           %The i-th sample sequence
                Ch = char(CHROMS(1,i));                         %The i-th sample sequence
                La = LABELS(1,i);                               %The i-th sample sequence
                L  = length(Se);                                %The length of this sequence
                if L>=81
                    Number_Length_81=L-81+1;                    %The number of fragments of length 81 in the i-th sample 
                    for j=1:Number_Length_81
                        fragment_j=Se(j:j+80);
                        NEW_SEQS.add(   fragment_j );
                        NEW_LABELS.add( La );
                        NEW_CHROMS.add( Ch );
                    end   
                    t=t+Number_Length_81;
                else
                    NEW_SEQS.add(   Se );
                    NEW_LABELS.add( La );
                    NEW_CHROMS.add( Ch );
                end
                waitbar(i/shuff_len, f , sprintf("%s",Se));
            end
            close(f)

            NEW_SEQS   = cell(toArray(NEW_SEQS))';
            NEW_LABELS = cell(toArray(NEW_LABELS))';
            NEW_CHROMS = cell(toArray(NEW_CHROMS))';

            SEQS   = strings(1, length(NEW_SEQS) );
            LABELS = zeros(1, length(NEW_LABELS) );
            CHROMS = strings(1, length(NEW_CHROMS) );

            for i=1:shuff_len 
               SEQS(1, i)   = char(NEW_SEQS(i));
               CHROMS(1, i) = char(NEW_CHROMS(i));
               LABELS(1, i) = cell2mat(NEW_LABELS(i));
            end

            disp(SEQS(1:5))

            fprintf('\nNEW SHAPES: \n\nCHROM SHAPE: (%s) \nSEQ SHAPE: (%s) \nLABELS: (%s) \n', num2str(size(CHROMS)), num2str(size(SEQS)), num2str(size(LABELS)));
            fprintf('\nSAMPLES: \n\nCHROM: %s \nSEQ: %s \nLABEL: %d \n', CHROMS(1), SEQS(1), LABELS(1) ); 

            SHUFFLE_T  = table( CHROMS.', SEQS.', LABELS.', 'VariableNames',{'CHROM','SEQ', 'LABEL'} );
            disp( head(SHUFFLE_T) );
            fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);

            disp("___________________SAVING TO DISK___________________");
            save(  bacteria_output_file_SEQS, 'SEQS'); 
            save(  bacteria_output_file_LABELS, 'LABELS');
            save(  bacteria_output_file_CHROMS, 'CHROMS');  
            save(  bacteria_output_file_SUMMARY, 'SHUFFLE_T');  
        else    
            disp("___________________LOADING FROM DISK___________________");
            load(  bacteria_output_file_SEQS    ); 
            load(  bacteria_output_file_LABELS  );
            load(  bacteria_output_file_CHROMS  );  
            load(  bacteria_output_file_SUMMARY );
            
            fprintf('\nSHAPES: \n\nCHROM SHAPE: (%s) \nSEQ SHAPE: (%s) \nLABELS: (%s) \n', num2str(size(CHROMS)), num2str(size(SEQS)), num2str(size(LABELS)));
            fprintf('\nSAMPLES: \n\nCHROM: %s \nSEQ: %s \nLABEL: %d \n', CHROMS(1), SEQS(1), LABELS(1) ); 
            
            if( preprocess_data == true )
                disp("___________________PRE-PROCESSING___________________");
                if(generate_KNN_data == true)
                    tic
                    fprintf('\n(1\5) PROCESSING KNN\n');
                    PPT1 = KNN(CHROMS,SEQS);
                    save( strcat(bacteria_output_dir, "/KNN.mat"), "PPT1");
                    fprintf('\n\nTIME ELAPSED (min): %.3f\n\n', toc/60);
                else
                    fprintf('\n(1\5) LOAD KNN\n');
                    load(strcat(bacteria_output_dir, "/KNN.mat"))
                end
                
                tic
                fprintf('(2\5) PROCESSING BPB');
                PPT2 = BPB(CHROMS,SEQS);
                fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);
                tic
                fprintf('(3\5) PROCESSING DNC');
                PPT3 = DNC(CHROMS,SEQS);
                fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);
                tic
                fprintf('(4\5) PROCESSING MNC');
                PPT4 = MNC(CHROMS,SEQS);
                fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);
                tic
                fprintf('(5\5) PROCESSING DAC');
                PPT5 = DAC(CHROMS,SEQS);
                fprintf('\n\nTIME ELAPSED (seconds): %.3f\n\n', toc);
    
                X = [PPT1,PPT2,PPT3,PPT4,PPT5]; 
                Y = LABELS';
    
                fprintf('\nX: (%s) \nY: (%s) \n', num2str(size(X)), num2str(size(Y)));
    
    
    
                fprintf("___________________SAVING DATA___________________\n\nX: %s\nY: %s", bacteria_output_file_X, bacteria_output_file_Y );
                save(  bacteria_output_file_X, 'X');
                save(  bacteria_output_file_Y, 'Y');
            end
        end
    end 
end
