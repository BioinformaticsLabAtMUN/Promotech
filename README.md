# Promotech

Description

## Requirements
Description
1. 
2.
3.

## Commands

1. `-V, --version` - Print Promotech's latest version.
2. `-GUI, --gui`   - Use the interactive GUI interface for predicting 40 nucleotide sequences. The interface do not work for whole genome prediction.
3. `-F, --fasta`   - Specify the location of the sequences or whole genome FASTA file location in disk. This command is used together with the `-P, --predict-sequences` or `-PG, --parse-genome` arguments.
4. `-M, --model` - Indicates the type of model used for prediction and the target output data type used during the genome parsing stage. The default value is `RF-HOT`.
   - The available options for **40nt sequences prediction** are `RF-HOT`, `RF-TETRA`, `GRU`, `LSTM`. 
   - The available options for **whole genome parsing and prediction** are `RF-HOT`, `GRU`, `LSTM`. 
5. `-TS, --test-samples` - Used for testing purposes during the genome parsing stage. A whole genome can be made of 4 million+ nucleotides and can take hours, depending of your system configuration to parse and predict. This command limits the number of sequences the sliding window cuts from the genome. It is used only with the `-PG, --parse-genome` argument.
5. `-PG, --parse-genome` - Use a sliding window to cut 40 nucleotide sequences from the whole genome in forward and reverse strand. The files are then saved to the "results" folder with a ".data" format and with the name of the model type, i.e. "RF-HOT.data" and "RF-HOT-INV.data". 
   - The **mandatory** argument used with this command is `-F, --fasta`. 
   - The **optional** arguments used with this command are `-M, --model`, and `--TS, --test-samples`. 
6. `-G, --predict-genome` -  This command uses the files generated using the `-PG, --parse-genome` argument and located in the "results" folder. 
   - The **optional** argument used with this command is `-M, --model`. Make sure to match the same model type used during the parsing stage.
7. `-P, --predict-sequences` - 
   - The **mandatory** argument used with this command is `-F, --fasta`. 
   - Make sure that the FASTA file has only 40-nt sequences as shown in the example below. If you require to use longer sequences, use the `-PG, --parse-genome` and `-G, --predict-genome` commands.
   ```
    >chrom1
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    >chrom2
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    >chrom3
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   ```



## Predicting Promoters Examples
### 40 Nucleotide Sequences
Description

### Whole Genome
Description


1. Parse the genome FASTA File

`python promotech.py -PG -F tests/genome/ECOLI_2.fasta` 

or smaller number of sequences can be used for testing purposes using the **--test-samples, -TS** parameter.

`python promotech.py -PG -TS 20000 -F tests/genome/ECOLI_2.fasta` 

2. Predict Promoter Sequences Using the Parsed Sequences

`python promotech.py -G `
