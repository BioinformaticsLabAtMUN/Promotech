# Promotech

Description

## Requirements
Description
1. Download and Install Anaconda or Miniconda from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). 
   - `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh`
   - `bash miniconda.sh`
2. Install conda environment from the prebuilt environment YAML file. 
   - `conda env create -f promotech_env.yml`
3. Activate environment
   - `conda activate promotech_env`

## Commands

1. `-v, --version` - Print Promotech's latest version.
2. `-gui, --gui`   - Use the interactive GUI interface for predicting 40 nucleotide sequences. The interface do not work for whole genome prediction.
3. `-f, --fasta`   - Specify the location of the sequences or whole genome FASTA file location in disk. This command is used together with the `-s, --predict-sequences` or `-PG, --parse-genome` arguments.
4. `-m, --model` - Indicates the type of model used for prediction and the target output data type used during the genome parsing stage. The default value is `RF-HOT`.
   - The available options for **40nt sequences prediction** are `RF-HOT`, `RF-TETRA`, `GRU`, `LSTM`. 
   - The available options for **whole genome parsing and prediction** are `RF-HOT`, `GRU`, `LSTM`. 
5. `-ts, --test-samples` - Used for testing purposes during the genome parsing stage. A whole genome can be made of 4 million+ nucleotides and can take hours, depending of your system configuration to parse and predict. This command limits the number of sequences the sliding window cuts from the genome. It is used only with the `-pg, --parse-genome` argument.
5. `-pg, --parse-genome` - Use a sliding window to cut 40 nucleotide sequences from the whole genome in forward and reverse strand. The files are then saved to the "results" folder with a ".data" format and with the name of the model type, i.e. "RF-HOT.data" and "RF-HOT-INV.data". 
   - The **mandatory** argument used with this command is `-f, --fasta`. 
   - The **optional** arguments used with this command are `-m, --model`, and `--ts, --test-samples`. 
6. `-g, --predict-genome` -  This command uses the files generated using the `-pg, --parse-genome` argument and located in the "results" folder. 
   - The **optional** argument used with this command is `-m, --model`. Make sure to match the same model type used during the parsing stage.
7. `-s, --predict-sequences` - Used to parse and predict 40nt sequences from a FASTA file.
   - The **mandatory** argument used with this command is `-f, --fasta`. 
   - Make sure that the FASTA file has only 40-nt sequences as shown in the example below. If you require to use longer sequences, use the `-pg, --parse-genome` and `-g, --predict-genome` commands.
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


1. Parse the whole genome FASTA File

`python promotech.py -PG -F examples/genome/ECOLI_2.fasta` 

or smaller number of sequences can be used for testing purposes using the **--test-samples, -TS** parameter.

`python promotech.py -PG -TS 20000 -F examples/genome/ECOLI_2.fasta` 

2. Predict promoter sequences using the parsed sequences

`python promotech.py -G `
