import pandas as pd
import numpy  as np 
import progressbar
import traceback

def charToBinary(argument):
  I = np.identity(4)
  switcher = {
      "A" : I[0],
      "G"	: I[1],
      "C"	: I[2],
      "T"	: I[3],
  }
  return switcher.get(argument)
def stringToBinary(string):
  tokenArray = list(string)
  result = list()
  for s in tokenArray:
    result = np.append( result , charToBinary(s))
  return result
def getHotFeatures(seqs, shape):
  print("GENERATING FEATURES FOR LABEL ", seqs.shape )
  data = np.array([]).reshape(shape)
  bar  = progressbar.ProgressBar(max_value=len( seqs ))
  for i , seq in enumerate( seqs ):
    binarySeq         = stringToBinary(seq)
    if(binarySeq.shape[0] < 160 ):
      continue
    data = np.vstack( [ data, binarySeq ] )
    bar.update(i)
  print( data.shape )
  return data
def fastaToHotEncoding(seqs):
  I = np.identity(4)
  columns = list([char for char in ("AGCT"*40) ]) 
  data = getHotFeatures(seqs, shape=(0, 160))
  return pd.DataFrame(data, columns=columns, dtype='int32')

def fastaToTetraNucletideDic(seqArray, label):
    pass

def parseSequences(seqs, sequence_length=40, slop_mode="middle", tokenizer_path=None, data_type="RF-HOT" ):
    if( len(seqs) == 0 ):
      print("No sequences available.", seqs)
      return
    sample, sample_len     = seqs[0], len(seqs[0])
    seqs = np.array([ s.upper() for s in seqs ])
    print("INPUT SHAPE {} \nSAMPLE WITH LEN {}: \n{}".format( len(seqs), sample_len, sample))
    if(sample_len > sequence_length):
        seqs = np.array([  s[0:sequence_length] for s in seqs ])
        print("OUTPUT SAMPLE WITH LEN {}: \n{}".format( seqs.shape, seqs[0] ))
    if data_type == "RF-HOT":
        data_df = fastaToHotEncoding( seqs )
    elif data_type == "RF-TETRA":
        data_df = fastaToTetraNucletideDic( seqs, None )
    elif data_type == "RNN":
        if( tokenizer_path  != None ):
          tokenizer          = pickle.load(open( tokenizer_path, 'rb' ))
        #   tetra_seq_list     = fasta_to_tetranucleotide_list(seq_df.values, np.array([]) )
        #   tetra_seq_list_str = tetranucleotide_list_to_string_list( tetra_seq_list.iloc[:, :-1].values )
        #   data_df            = pd.DataFrame( tokenizer.texts_to_sequences(tetra_seq_list_str) )
    print(data_df.head())
    return data_df
