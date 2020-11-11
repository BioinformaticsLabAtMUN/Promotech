___________________________________________________________________
bTSSfinder README

    INTRODUCTION
    PREREQUISITES
    HOW TO RUN
    EXAMPLE
    OPTIONS EXPLAINED
    BED FORMAT
___________________________________________________________________

AUTHORS:	I.A.Shahmuradov, Rozaimi Mohamad Razali, Salim Bougouffa, Aleksandar Radovanovic and Vladimir Bajic
LAST UPDATE:    08 June 2016
VERSION:        1.2016
ACCESS:         ..........................

___________________________________________________________________
INTRODUCTION:
___________________________________________________________________

bTSSfinder was developed in the Computational Bioscience Research Center (CBRC)
     at the King Abdullah University of Science and Technology (KAUST), Kingdom of Saudi Arabia.
     It aims to predict the transcription start sites (TSSs) of 5 sigma factor classes
     for Escherichia coli (sigma 70, 38, 32,28 and 24) and Cyanonacteria (sigma A, C, H, F and G).

Input query file: a single or multiple sequences in FASTA format.
     Maximum allowed length for each sequence: 25,000,000 bp.

Output file is in three formats: (1) classic Text format, (2) GFF3 format, (3) BED format.

___________________________________________________________________
PREREQUISITES:
___________________________________________________________________

bTSSfinder can run on Linux or Unix (MacOS).
    - GFORTRAN "gfortran" is required.

___________________________________________________________________
HOW TO RUN:
___________________________________________________________________

*** REQUIRED ***
   You MUST set up the environmental variable bTSSfinder_Data,
   Data is included in the program

   For example:
        setenv bTSSfinder_Data /path/to/installation/Data
   or
        export bTSSfinder_Data=/path/to/installation/Data
   or
	     Add the following to your ~/.bashrc or ~/.bash_profile

      	export PATH=/path/to/installation:$PATH
      	export bTSSfinder_Data=/path/to/installation/Data/:$bTSSfinder_Data

___________________________________________________________________
EXAMPLE:
___________________________________________________________________

	/path/to/installation/bTSSfinder -i example_Cyanob.fasta -o example_Cyanob -h 2 -r 0 -t c
___________________________________________________________________
OPTIONS EXPLAINED
___________________________________________________________________

 bTSSfinder -i <parI> - Input fasta file  * REQUIRED
      [-p <parP> - print (y or Y) or not print (n or N) the query sequence(s): Default: n]
      [-o <parO> - Prefix for output file; Default: bTSSfinder]
      [-x <parX> - TSS search window; >=50,<=300; Default: 300]
      [-n <parN> - Comma-separated A,C,G,T frequencies in the reference genome; Default: 0.25,0.25,0.25,0.25]
      [-a <parA> - Neural Network Threshold for sigma70 or sigmaA promoters;  -2<=parA<=2, Default: preset in training ]
      [-b <parB> - Neural Network Threshold for sigma38 or sigmaC promoters;  -2<=parB<=2, Default: preset in training ]
      [-d <parD> - Neural Network Threshold for sigma32 or sigmaH promoters;  -2<=parD<=2, Default: preset in training ]
      [-e <parE> - Neural Network Threshold for sigma28 or sigmaF promoters;  -2<=parE<=2, Default: preset in training ]
      [-f <parF> - Neural Network Threshold for sigma24 or sigmaG promoters;  -2<=parD<=2, Default: preset in training ]
      [-h <parH> - 1 to search on the sense strand OR 2 for both strands, Default: 1 ]
      [-t <parT> - Taxon (e or E - for E. coli, and c or C - for Cyanobacteria, Default: E ]
      [-c <parC> - Search Criterion; if parC =
                                     70 ... search for sigma70 promoters (E. coli)
                                     38 ... search for sigma38 promoters (E. coli)
                                     32 ... search for sigma32 promoters (E. coli)
                                     28 ... search for sigma28 promoters (E. coli)
                                     24 ... search for sigma24 promoters (E. coli)
                                      A ... search for sigmaA promoters (cyanobacteria)
                                      C ... search for sigmaC promoters (cyanobacteria)
                                      F ... search for sigmaF promoters (cyanobacteria)
                                      G ... search for sigmaG promoters (cyanobacteria)
                                      H ... search for sigmaH promoters (cyanobacteria)
                                      1 ... search for the highest ranking promoter regardless of class in the chosen Taxon (option -t)
                                      5 ... [Default], report all promoter classes in the chosen Taxon (-t)

___________________________________________________________________
BED FORMAT:
___________________________________________________________________

The output file in the BED format has the following fields:

(1) The DNA segment identifier usually the fasta ID (e.g. chr1, gene1).
(2) The start position of the predicted TSS.
(3) The end position of the predicted TSS.
(4) Feature name - TSS#.
(5) Score - 0.
(6) Strand '+' or '-'.
(7) thickStart - TSS# start position.
(8) thickEnd TSS# End position, the same as (7).
(9) An RGB of the form R,G,B (set to 255,0,0).
(10) The number of main promoter blocks (in our case =2, -35 and -10 boxes).
(11) A comma-separated list of block sizes(in our case = 6,6).
(12) A comma-separated list of block start positions (calculated relative to the starting position of the query).
