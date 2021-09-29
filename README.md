# PROMOTECH: A UNIVERSAL TOOL FOR PROMOTER DETECTION IN BACTERIAL GENOMES

Promotech is a machine-learning-based classifier trained to generate a model that generalizes and detects promoters in a wide range of bacterial species.  During the study, two model architectures  were  tested,  Random  Forest  and  Recurrent  Networks. The  Random Forest model, trained with promoter sequences with a binary encoded representation of  each  nucleotide,  achieved  the  highest  performance  across  nine  different  bacteria and was able to work with short 40bp sequences and whole bacterial genomes using a sliding window.  The selected model was evaluated on a validation set of four bacteria not used during training.

## Requirements

1. Download and Install Anaconda or Miniconda from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). 
   - `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh` for *Ubuntu 20.04*
   - `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh` for *Mac OS Big Sur V11.3*
   - `bash miniconda.sh`
2. Install conda environment from the prebuilt environment YAML file. 
   - `conda env create -f promotech_env.yml` for *Ubuntu 20.04*
   - `conda env create -f promotech_mac_env.yml` for *Mac OS Big Sur V11.3*
   - **Note:** The environment was made on Ubuntu 20.04, different versions of the packages could be required for different operating systems.
3. Activate environment
   - `conda activate promotech_env`
4. Download the models from [here](http://www.cs.mun.ca/~lourdes/public/PromoTech_models/)
5. Uncompress the two Random Forest models with *.zip* format and save them at the [*models* folder](models/). The resulting files are *models/RF-HOT.model* and *models/RF-TETRA.model*
6. Make sure all files required to run Promotech are available in promotech's folder

**Note:** A minimum of 24 GB of RAM memory is recommended to run the RF-HOT, LSTM, and GRU model on a whole-genome. Parsing a whole-genome to the RF-TETRA model's input format can produce the python "Memory Error" due to the high complexity and high RAM-memory demand required to obtain the tetra-nucleotide frequencies for millions of sequences in forward and inverse strand. An example of this process is shown in the examples section below. All models can run on lower-end systems, with at least 8GB of RAM, when predicting FASTA files with hundreds or thousands of sequences, 40 nt in length. 

The examples in the section below were tested in a desktop computer with the following specifications:

- **Processor      :** Intel(R) Core(TM) i5-9300H CPU @ 2.40GHz 2.40 GHz 
- **RAM            :** 24.0 GB (23.8 GB usable)
- **System Type    :** 64-bit Ubuntu 20.04 LTS
- **Graphic Memory :** NVIDIA GeForce RTX 2060 6GB GDDR6
- **Python Version :** Python 3.6


## Commands

1. `-v, --version` - Print Promotech's latest version.
2. `-gui, --gui`   - Use the interactive GUI interface for predicting 40 nucleotide sequences. The interface does not work for whole-genome prediction.
3. `-f, --fasta`   - Specify the location of the sequences or whole-genome FASTA file location in disk. This command is used together with the `-s, --predict-sequences` or `-PG, --parse-genome` arguments.
4. `-m, --model` - Indicates the type of model used for prediction and the target output data type used during the genome parsing stage. The default value is `RF-HOT`.
   - The available options for **40nt sequences prediction** are `RF-HOT`, `RF-TETRA`, `GRU`, `LSTM`. 
   - The available options for **whole-genome parsing and prediction** are `RF-HOT`, `GRU`, `LSTM`. 
5. `-ts, --test-samples` - Used for testing purposes during the genome parsing stage. A whole-genome can be made of 4 million+ nucleotides and can take hours, depending on your system configuration to parse and predict. This command limits the number of sequences the sliding window cuts from the genome. It is used only with the `-pg, --parse-genome` argument.
5. `-pg, --parse-genome` - Use a sliding window to cut 40 nucleotide sequences from the whole-genome in forward and reverse strand. The files are then saved to the "results" folder with a "[MODEL-TYPE].data" format, where MODEL-TYPE is the name of the model's desired input format, i.e. "RF-HOT.data" and "RF-HOT-INV.data". 
   - The **mandatory** argument used with this command is `-f, --fasta`. 
   - The **optional** arguments used with this command are `-m, --model`, and `--ts, --test-samples`. 
6. `-g, --predict-genome` -  This command uses the files generated using the `-pg, --parse-genome` argument and located in the "results" folder. 
   - The **optional** argument used with this command is `-m, --model`. Make sure to match the same model type used during the parsing stage.
7. `-s, --predict-sequences` - Used to parse and predict 40nt sequences from a FASTA file.
   - The **mandatory** argument used with this command is `-f, --fasta`. 
   - Make sure that the FASTA file has only 40-nt sequences as shown in the example below. If you require to use longer sequences, use the `-pg, --parse-genome` and `-g, --predict-genome` commands.
  <br />
  
  ```
    >chrom1
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    >chrom2
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    >chrom3
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   ```



## Examples

### 40 Nucleotide Sequences
1. Parse and predict using the `--predict-sequence, -s` command and specify the FASTA file using `--fasta, -f`. The FASTA file should only include 40nt sequences. If you require to predict longer sequences, use the "whole-genome" commands. An example of the command-line output obtained when running this command is found [HERE](examples/sequences/example.log.txt)

`python promotech.py -s -m RF-HOT -f examples/sequences/test.fasta -o results`

### Whole-Genome

1. Parse the whole-genome in the FASTA file (a single sequence is expected) by using `--parse-genome, -pg` and specifying the file using `--fasta, -f` . A smaller subset of the sliding window sequences can be used for testing purposes using the `--test-samples, -ts` parameter.  An example of the command-line output obtained when running this command is found [HERE](examples/genome/parse.example.log.txt)

`python promotech.py -pg -m RF-HOT -f examples/genome/ECOLI_2.fasta -o results` 

or 

`python promotech.py -pg -ts 20000 -m RF-HOT -f examples/genome/ECOLI_2.fasta -o results` 

- **Note:** Running one of the following commands will use a sliding window of 40nt size and 1nt step, pre-processed the sequences to meet the specified model's input requirement and create two files, **results/[MODEL_TYPE].data** and **results/[MODEL_TYPE]-INV.data**, for forward and inverse strand, where MODEL_TYPE specifies the type of model that will later be used for assessing the pre-processed 40nt sequences. **Do not delete the 'results' folder or the '.data' files**, because they will be used in the next step.
- **Note:** For comparison, the pipeline configured to generate data for the RF-HOT model, took 35 minutes and 42 seconds to cut 4,639,634 forward and 4,639,634 inverse sequences from the *E. coli* (NC_000913.2) genome with 4,639,675 nucleotides in length, pre-processed them to hot-encoded binary format, save them to two binary files and each file was 5.8 GB in size. During this time, it maintained around 18.5/24GB of RAM exclusively for the python running process.

2. Predict promoter sequences using the parsed sequences using `--predict-genome, -g`, assign a threshold using `--threshold, -t`, and select a model using `--model, -m`. The default threshold, and model are 0.5, and RF-HOT, respectively. An example of the command-line output obtained when running this command is found [HERE](examples/genome/predict.example.log.txt)

`python promotech.py -g -t 0.6 -i results -o results`

- **Note:** This command expects the user to have used the `--parse-genome, -pg` command before to generate the pre-processed sequences from the bacterial genome and stored in the files **results/[MODEL_TYPE].data** and **results/[MODEL_TYPE]-INV.data**. 
- **Note:** For comparison, it took 1 hour, 5 minutes, and 27 seconds to predict both, forward and inverse strand batches, each with 4,639,634 pre-processed sequences, with a total of 9,279,268 sequences as input and an output of 55,002 promoters sequences with a score above the 0.5 threshold.

## Cite

If you use Promotech please cite:

Promotech: A general tool for bacterial promoter recognition. Ruben Chevez-Guardado and Lourdes Peña-Castillo. [bioRxiv 2021.07.16.452684](https://doi.org/10.1101/2021.07.16.452684)
