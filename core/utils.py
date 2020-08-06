import pandas as pd
import numpy  as np 
from sklearn.model_selection import train_test_split
from skbio import Sequence
from itertools import product
import progressbar
import traceback

def readData(positives_filename, negatives_filename):
  p_df = pd.read_csv(positives_filename, sep='\t', header=None)
  n_df  = pd.read_csv(negatives_filename, sep='\t', header=None)
  p_df = p_df[p_df[0].str.contains(">")==False]
  n_df = n_df[n_df[0].str.contains(">")==False]
  return p_df, n_df

def all_tetra_subsets(ss=["A", "T", "G", "C"]):
  subsets = list(['AAAA', 'AAAT', 'AAAG', 'AAAC', 'AATA', 'AATT', 'AATG', 'AATC', 'AAGA', 'AAGT', 'AAGG', 'AAGC', 'AACA', 'AACT', 'AACG', 'AACC', 'ATAA', 'ATAT', 'ATAG', 'ATAC', 'ATTA', 'ATTT', 'ATTG', 'ATTC', 'ATGA', 'ATGT', 'ATGG', 'ATGC', 'ATCA', 'ATCT', 'ATCG', 'ATCC', 'AGAA', 'AGAT', 'AGAG', 'AGAC', 'AGTA', 'AGTT', 'AGTG', 'AGTC', 'AGGA', 'AGGT', 'AGGG', 'AGGC', 'AGCA', 'AGCT', 'AGCG', 'AGCC', 'ACAA', 'ACAT', 'ACAG', 'ACAC', 'ACTA', 'ACTT', 'ACTG', 'ACTC', 'ACGA', 'ACGT', 'ACGG', 'ACGC', 'ACCA', 'ACCT', 'ACCG', 'ACCC', 'TAAA', 'TAAT', 'TAAG', 'TAAC', 'TATA', 'TATT', 'TATG', 'TATC', 'TAGA', 'TAGT', 'TAGG', 'TAGC', 'TACA', 'TACT', 'TACG', 'TACC', 'TTAA', 'TTAT', 'TTAG', 'TTAC', 'TTTA', 'TTTT', 'TTTG', 'TTTC', 'TTGA', 'TTGT', 'TTGG', 'TTGC', 'TTCA', 'TTCT', 'TTCG', 'TTCC', 'TGAA', 'TGAT', 'TGAG', 'TGAC', 'TGTA', 'TGTT', 'TGTG', 'TGTC', 'TGGA', 'TGGT', 'TGGG', 'TGGC', 'TGCA', 'TGCT', 'TGCG', 'TGCC', 'TCAA', 'TCAT', 'TCAG', 'TCAC', 'TCTA', 'TCTT', 'TCTG', 'TCTC', 'TCGA', 'TCGT', 'TCGG', 'TCGC', 'TCCA', 'TCCT', 'TCCG', 'TCCC', 'GAAA', 'GAAT', 'GAAG', 'GAAC', 'GATA', 'GATT', 'GATG', 'GATC', 'GAGA', 'GAGT', 'GAGG', 'GAGC', 'GACA', 'GACT', 'GACG', 'GACC', 'GTAA', 'GTAT', 'GTAG', 'GTAC', 'GTTA', 'GTTT', 'GTTG', 'GTTC', 'GTGA', 'GTGT', 'GTGG', 'GTGC', 'GTCA', 'GTCT', 'GTCG', 'GTCC', 'GGAA', 'GGAT', 'GGAG', 'GGAC', 'GGTA', 'GGTT', 'GGTG', 'GGTC', 'GGGA', 'GGGT', 'GGGG', 'GGGC', 'GGCA', 'GGCT', 'GGCG', 'GGCC', 'GCAA', 'GCAT', 'GCAG', 'GCAC', 'GCTA', 'GCTT', 'GCTG', 'GCTC', 'GCGA', 'GCGT', 'GCGG', 'GCGC', 'GCCA', 'GCCT', 'GCCG', 'GCCC', 'CAAA', 'CAAT', 'CAAG', 'CAAC', 'CATA', 'CATT', 'CATG', 'CATC', 'CAGA', 'CAGT', 'CAGG', 'CAGC', 'CACA', 'CACT', 'CACG', 'CACC', 'CTAA', 'CTAT', 'CTAG', 'CTAC', 'CTTA', 'CTTT', 'CTTG', 'CTTC', 'CTGA', 'CTGT', 'CTGG', 'CTGC', 'CTCA', 'CTCT', 'CTCG', 'CTCC', 'CGAA', 'CGAT', 'CGAG', 'CGAC', 'CGTA', 'CGTT', 'CGTG', 'CGTC', 'CGGA', 'CGGT', 'CGGG', 'CGGC', 'CGCA', 'CGCT', 'CGCG', 'CGCC', 'CCAA', 'CCAT', 'CCAG', 'CCAC', 'CCTA', 'CCTT', 'CCTG', 'CCTC', 'CCGA', 'CCGT', 'CCGG', 'CCGC', 'CCCA', 'CCCT', 'CCCG', 'CCCC'])   
  return subsets

def tetraToHotEncoding(seq_arr):
  I = np.identity(4)
  columns = list([char for char in ("AGCT"*40) ]) 
  seq_arr = [ seq.strip() for seq in seq_arr ]
  print("\n\n ", "_"*10 , "TETRA-NUCLEOTIDE SEQ TO HOT ENCODING CONVERSION" , "_"*10 , "\n\n")

  print("TetraToHotEncoding Input Sample: ", seq_arr[0] )
  print("TetraToHotEncoding Cols: "," A: ", I[0], " G: ", I[1]," C: ",I[2]," T: ",I[3] )

  hot_seq_arr = np.empty( (len(seq_arr), 160) , dtype=int )
  bar         = progressbar.ProgressBar(max_value=len(seq_arr))
  for i, seq in enumerate(seq_arr):
    c_c = 0
    for bp in list(seq):
      chars_data = charToBinary(bp) 
      hot_seq_arr[i, c_c+0] = chars_data[0]
      hot_seq_arr[i, c_c+1] = chars_data[1]
      hot_seq_arr[i, c_c+2] = chars_data[2]
      hot_seq_arr[i, c_c+3] = chars_data[3]
      c_c += 4
    bar.update(i)
  return hot_seq_arr

def promoterToTetraFreq(seq_arr):
  print("\n\n ", "_"*10 , "PROMOTER SEQ TO TETRA FRQUENCIES CONVERSION" , "_"*10 , "\n\n")
  print("Input Sample: ", seq_arr[0] )
  tetra_n        = np.empty( len(seq_arr), dtype=object )
  nucleotide_counter = 0
  nucleotide_list = ["U","R","Y","N","W","S","M","K","B","H","D","V" ]
  columns=all_tetra_subsets()
  bar         = progressbar.ProgressBar(max_value=len(seq_arr))
  for i, s in enumerate(seq_arr):
    if any( nucleotide in s for nucleotide in nucleotide_list ):
      nucleotide_counter = nucleotide_counter + 1
      continue
    if(len(s) < 40 ):
      nucleotide_counter = nucleotide_counter + 1
      continue
    seq_dict = Sequence( s ).kmer_frequencies(4, overlap=True, relative =True)
    tetra_n[i] = seq_dict
    bar.update(i)
  print("nucleotide discarded counter: ",nucleotide_counter , "\n\n")
  tetra_n_df = pd.DataFrame.from_records(tetra_n, columns=columns ) 
  print("Checking for missing columns")
  col_added = 0
  for col in columns: 
    if col not in tetra_n_df.columns:
      tetra_n_df[col] = 0
      col_added+=1
  print("Reordering columns and looking for nan values. Columns added: ", col_added)
  tetra_n_df = tetra_n_df.replace(np.nan, 0)
  print("\n\nDONE. Tetra Freq Shape & Sample: ",tetra_n_df.shape , "\n", tetra_n_df.head() , "\n")
  return tetra_n_df
  
def fastaToTetraNucletideDic(seqArray, label):
  print("FASTA TO TETRA-NUCLEOTIDE FOR LABEL: ", label, " - ", seqArray.shape )
  tetra_n = list()
  nucleotide_counter = 0
  #Only "A","G","C","T" combinations accepted
  nucleotide_list = ["U","R","Y","N","W","S","M","K","B","H","D","V" ]
  sequence_arr = np.empty(len(seqArray), dtype=object)
  bar  = progressbar.ProgressBar(max_value=len(seqArray))
  for i, s in enumerate(seqArray):
    if any( nucleotide in s[0] for nucleotide in nucleotide_list ):
      nucleotide_counter = nucleotide_counter + 1
      sequence_arr[i] = dict()
      continue
    if(len(s[0]) < 40 ):
      nucleotide_counter = nucleotide_counter + 1
      sequence_arr[i] = dict()
      continue
    seq = Sequence(s[0]).kmer_frequencies(
        4, overlap=True, relative =True)
    sequence_arr[i] = seq
    bar.update(i)
  tetra_n_df = pd.DataFrame.from_records(sequence_arr)
  tetra_n_df = tetra_n_df[all_tetra_subsets()]
  tetra_n_df = tetra_n_df.replace(np.nan, 0)
  labels_df = pd.Series([ label ] *  len(tetra_n_df.index)  )
  data = pd.concat([tetra_n_df, labels_df.rename("label")], axis=1)   
  print("nucleotide discarded counter: ",nucleotide_counter, "\n\n")
  print("TETRA FREQ DATA SHAPE: ", data.shape, "\n\n")
  print(data.head())
  return data 

def fastaToCharArray(seqArray, label):
  charArray = [list(s[0]) for s in seqArray]
  for c in range(len(charArray)):
    charArray[c].append(label)
  return charArray

def joinPositiveAndNegative(positive_tetra_n, negative_tetra_n ):
  print("P: ", positive_tetra_n.shape)
  print("N: ", negative_tetra_n.shape)
  join_tetra_n = pd.concat([ positive_tetra_n,  negative_tetra_n ],  ignore_index=True )
  return join_tetra_n

def generateTrainAndTestSplit(dataset):
  X = dataset[:, :-1]
  y = dataset[:,  -1]
  print("generateTrainAndTestSplit Y", y)
  print("Complete Data+Label Shape: ", dataset.shape, " Features shape ", X.shape," Labels Shape: ",  y.shape, " Pos: ", sum(map(lambda x : x == 1, y))  , "Neg: ", sum(map(lambda x : x == 0, y))   )
  return train_test_split( 
      X, 
      y, 
      test_size=0.10, 
      random_state=42,
      shuffle= True,
      stratify=y
  )

################################## HOT-ENCODING FUNCTIONS ##############################################
def charToBinary(argument):
  #https://www.genome.jp/kegg/catalog/codes1.html
  # I = np.identity(16)
  I = np.identity(4)
  switcher = {
      "A" : I[0],
      "G"	: I[1],
      "C"	: I[2],
      "T"	: I[3],
      # "U"	: I[4],
      # "R"	: I[5],
      # "Y"	: I[6],
      # "N"	: I[7],
      # "W"	: I[8],
      # "S"	: I[9],
      # "M"	: I[10],
      # "K" : I[11],
      # "B" : I[12],
      # "H"	: I[13],
      # "D"	: I[14],
      # "V" : I[15]
  }
  return switcher.get(argument)

def sequenceToBinary( sequence ):
  return np.array([ charToBinary(nucleotide) for nucleotide in sequence ]).flatten() #map( charToBinary, sequence ) 

def getHotFeatures(sequences):
  n_seqs   = len( sequences )
  if n_seqs <= 0:
    return
  sample   = sequences[0]
  data     = np.empty( (n_seqs, len(sample)*4) , dtype=np.ndarray)
  bar      = progressbar.ProgressBar(max_value=n_seqs)
  for i , seq in enumerate( sequences ):
    data[i, :] = sequenceToBinary(seq); bar.update(i)
  return data
def fastaToHotEncodingSequences(seqs, nucleotide_order="AGCT"):
  I = np.identity(4)
  columns = list([char for char in (nucleotide_order*40) ]) 
  data = getHotFeatures(seqs)
  return pd.DataFrame(data, columns=columns, dtype='int32')

#########################################################################################################
    
def stringToBinaryArray(string):
  tokenArray = list(string)
  result = list()
  for s in tokenArray:
    result = np.append( result , charToBinary(s))
  return result

def getFeaturesAndLabels(dataframe, label, shape):
  print("GENERATING FEATURES FOR LABEL ", label, " - ", dataframe.shape )
  data = np.array([]).reshape(shape)
  bar  = progressbar.ProgressBar(max_value=len( dataframe.values ))
  for i , seq in enumerate( dataframe.values ):
    binarySeq         = stringToBinaryArray(seq[0])
    binarySeqAndLabel = np.append(  binarySeq , label)
    # print(  "Append", data.shape, binarySeqAndLabel.shape , binarySeq)
    if(binarySeqAndLabel.shape[0] < 161 ):
      continue
    data = np.vstack( [ data, binarySeqAndLabel ] )
    bar.update(i)
  print( data.shape )
  return data
  
def fastaToHotEncoding(p_df, n_df):
  I = np.identity(4)
  columns = list([char for char in ("AGCT"*40) ]) 
  columns.append("label")
  p_data = getFeaturesAndLabels(p_df, label=1, shape=(0, 161))
  n_data = getFeaturesAndLabels(n_df, label=0, shape=(0, 161))
  return pd.DataFrame(p_data, columns=columns, dtype='int32'), pd.DataFrame(n_data, columns=columns, dtype='int32')

def fasta_to_tetranucleotide_list(p_data, n_data):
  template = all_tetra_subsets()
  columns=all_tetra_subsets()
  columns.append('label')
  data   = []
  print("\n fasta_to_tetranucleotide_list")
  if( n_data is not None and len(n_data) > 0 ):
    print("\n n_data", n_data.shape )
    print("\n n_data", n_data[0] )
  if( p_data is not None and len(p_data) > 0 ):
    print("\n p_data", p_data.shape )
    print("\n p_data", p_data[0] )
  datasets = [n_data, p_data]
  for key in [1, 0]:
    label = key
    for sequence_dic in datasets[key]:
      if(sequence_dic is None or len(sequence_dic) <= 0 ):
        continue
      tetra_neucleotide_list = list()
      for k in range(len(sequence_dic[0])):
        segment = sequence_dic[0][k :  k+4 if k+4 <= len(sequence_dic[0]) else len(sequence_dic[0]) ]
        if len(segment) == 4:
          tetra_neucleotide_list.append(segment)
      tetra_neucleotide_list.append(label)
      data.append(tetra_neucleotide_list)
  data = pd.DataFrame(data )
  return data

def tetranucleotide_list_to_string_list(data):
  tetra_list = np.empty(len(data), dtype=object )
  for i, tetras in enumerate(data): 
    tetra_list[i] = " ".join(tetras)
  return tetra_list

def testTry(fn, addErrMsg=""):
    try:
        fn()
    except Exception as e:
        print("Error:", addErrMsg,"-", repr(e))

def h1(logs_arr, msg):
  print( "\n\n", msg, "\n\n" )
  if(logs_arr != None ):
    logs_arr.append(msg)
  return logs_arr

def print_fn( msg, out_file ):
  print(msg)
  f = open(out_file, 'a')
  f.write(msg)
  f.close()
