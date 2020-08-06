import pandas as pd
import numpy  as np 
import os
import os.path
import traceback
import subprocess
import csv
from utils import readData, fastaToCharArray, fastaToHotEncoding, fastaToTetraNucletideDic
from utils import tetranucleotide_list_to_string_list, tetraToHotEncoding, promoterToTetraFreq
from utils import joinPositiveAndNegative, generateTrainAndTestSplit, h1
from utils import fasta_to_tetranucleotide_list, tetranucleotide_list_to_string_list
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.preprocessing.sequence import pad_sequences
import pickle
import progressbar
import warnings
import datetime
import time
from Bio.Seq import Seq 
import joblib
from pybedtools import BedTool
import pickle
warnings.filterwarnings("ignore")
from Bio import SeqIO

using_unbalanced = False

def createIndex():
  df = pd.DataFrame([
    #TRAINING
      ["ECOLI_2"      , "NC_000913.2"  , True],
      ["ECOLI"        , "NC_000913.3"  , True],
      ["HPYLORI"      , "NC_000915.1"  , True],
      ["HPYLORI_2"    , "NC_000915.1"  , True],
      ["C_JEJUNI"     , "NC_002163.1"  , True],
      ["CJEJUNI_2"    , "NC_002163.1"  , True],
      ["CJEJUNI_3"    , "NC_003912.7"  , True],
      ["CJEJUNI_4"    , "NC_009839.1"  , True],
      ["CJEJUNI_5"    , "NC_008787.1"  , True],
      ["SPYOGENE"     , "LR031521.1"   , True],   
      ["STYPHIRMURIUM", "NC_016810.1"  , True], 
      ["CPNEUMONIAE"  , "NC_000922.1"  , True],  
      ["SONEIDENSIS"  , "NC_004347.2"  , True],
      ["LINTERROGANS" , "NZ_LT962963.1", True],  
      ["SCOELICOLOR"  , "NC_003888.3"  , True],
      #VALIDATION
      ["MYCOBACTER"   , "NC_008596"   , False],
      ["CLOSTRIDIUM"  , "NC_010001.1" , False],
      ["RHODOBACTER_1", "NC_014034.1" , False],
      ["RHODOBACTER_2", "NC_014034.1" , False],
      ["BACILLUS"     , "CP002927.1" , False],
      # ["AERUGINOSA"   , "NC_002516.2" , False]
  ], columns=["BACTERIA", "ID", "IS_TRAINING"])
  return df

def saveIndex(df, dir):
  df.to_csv(dir)

def generateData(
  bacteria_index, save_csv=False, save_data=True, 
  out_dir="./data/promoters/", nextflow_path='./nextflow',
  nextflow_pipeline = "pipeline_unbalanced.nf", # 'pipeline_without_docker.nf'
  manually_balance_data = False
):
  if(using_unbalanced):
    print("GENERATING UNBALANCED DATA WITH RATIO 1:10")
  else:
    print("GENERATE DATA")

  # bacteriaDir = "./bacteria"
  bacteria_report = {}
  if bacteria_index is None:
    index = createIndex()
  else: 
    index = bacteria_index

  data_root = "./data/"
  if not os.path.exists(data_root):
      os.makedirs(data_root)

  w = csv.writer(open(data_root+"report.csv", "w"))
  vocab_size = None
  tokenizer = None

  hot_encoded_train_features = np.empty((0, 160), int)
  hot_encoded_train_labels   = np.empty((0,  ), int)
  hot_encoded_test_features = np.empty((0, 160), int)
  hot_encoded_test_labels   = np.empty((0,  ), int)
  hot_encoded_val_features = np.empty((0, 160), int)
  hot_encoded_val_labels   = np.empty((0,  ), int)


  tetra_freq_train_features = np.empty((0, 256), int)
  tetra_freq_train_labels   = np.empty((0,  ), int)
  tetra_freq_test_features = np.empty((0, 256), int)
  tetra_freq_test_labels   = np.empty((0,  ), int)
  tetra_freq_val_features = np.empty((0, 256), int)
  tetra_freq_val_labels   = np.empty((0,  ), int)


  rnn_token_train_features = np.empty((0, 37), int)
  rnn_token_train_labels   = np.empty((0,  ), int)
  rnn_token_test_features  = np.empty((0, 37), int)
  rnn_token_test_labels    = np.empty((0,  ), int)
  rnn_token_val_features  = np.empty((0, 37), int)
  rnn_token_val_labels    = np.empty((0,  ), int)
  global_rnn_complete       = np.empty((0, 37  ), int)

  start_time = datetime.datetime.now().time().strftime('%H:%M:%S')
  bar = progressbar.ProgressBar(max_value=len(index))
  for i, row in index.iterrows() :
    bacteria_start_time = datetime.datetime.now().time().strftime('%H:%M:%S')

    # print("\n\n", 20*"*", i+1, f". {row['BACTERIA']}", 20*"*" )
    print("\n\n {} {} {} {}".format(20*"*", i+1, row['BACTERIA'], 20*"*"  ) )
    #nextflow run main_pipeline.nf --bacteria ecoli && rsync outDir/ outDirOriginal/ -a --copy-links -v
    print("\n\n {} {} {} {}".format(20*"*", i+1, "NEXTFLOW DATA GENERATION", 20*"*"  ) )
    # print("\n\n", 10*"*", "NEXTFLOW DATA GENERATION",10*"*" )
    
    stderr = None      
    stdout = None
    
    if(nextflow_path is not None):
      print("\n\nGENERATING NEXTFLOW DATA USING PIPELINE: ", nextflow_pipeline, "\n\n")
      out = subprocess.Popen([
          nextflow_path,
          'run',
          nextflow_pipeline,   #'pipeline_without_docker.nf',   #    pipeline_unbalanced_without_docker.nf   'main_pipeline.nf',  
          '--bacteria',
          str(row['BACTERIA']),
      ],
      stdout=subprocess.PIPE,
      stderr=subprocess.STDOUT)
      stdout, stderr = out.communicate()
      error_msg = ""

      print("\n\nOUT: \n\n", stdout)
      print("\n\nERRORS: \n\n ",  stderr)

      bacteria_report[row['BACTERIA']] = {
        'stdout': stdout, 'stderr' : stderr
      }
    else:
      print("NEXTFLOW GENERATION SKIPPED.")
      


    if stderr == None:
      # print("\n\nConverting symlink to copy of file", row['BACTERIA'])
      # cmd = f"rsync outDir/{row['BACTERIA']} outDirOriginal/ -a --copy-links -v"
      
      if(nextflow_path is not None):
        cmd = "rsync "+out_dir+str(row['BACTERIA'])+" outDirOriginal/ -a --copy-links -v"
        out = os.popen(cmd).read()

      try:
        p_df, n_df = readData(
          out_dir+str(row['BACTERIA'])+"/positive.fasta", 
          out_dir+str(row['BACTERIA'])+"/negative.fasta") 
        
        p_bed_df = BedTool(out_dir+str(row['BACTERIA'])+"/positive.bed").to_dataframe()
        n_bed_df = BedTool(out_fir+str(row['BACTERIA'])+"/negative.bed").to_dataframe()
        p_bed_df["sequence"] = p_df.values
        n_bed_df["sequence"] = n_df.values
        p_bed_df["label"] = [1]*len(p_df)
        n_bed_df["label"] = [0]*len(n_df)
        dataset_df = pd.concat([p_bed_df, n_bed_df])
        print("SAVING DATASET: P {}  + N {} = {}".format( p_bed_df.shape, n_bed_df.shape, dataset_df.shape ))
        p_bed_df.to_csv(out_dir+str(row['BACTERIA'])+"/positive.csv")
        n_bed_df.to_csv(out_dir+str(row['BACTERIA'])+"/negative.csv")
        dataset_df.to_csv(out_dir+str(row['BACTERIA'])+"/dataset.csv")
        
        print("\n\n"+ 10*"*"+ "FASTA TO HOT ENCODING" + 10*"*" )
        print("P: {} N: {}".format(len(p_df) , len(n_df)))
        if(manually_balance_data and len(p_df) < len(n_df) ):
          print("Manually balancing Positives and Negatives. Decreasing Negatives from {} -> {}. Ratio {}:{}".format(
            len(n_df), len(p_df), 1, len(p_df) * 100 / len(n_df)
          ))
          n_df = n_df.sample(n=len(p_df))
          print("FINAL DATA SHAPES -> P: {} N : {}".format( p_df.shape, n_df.shape ))
          
        hot_p_data, hot_n_data = fastaToHotEncoding(p_df, n_df)
        hot_encoded_dataset_df = joinPositiveAndNegative(hot_p_data , hot_n_data )
        print("\n\n", hot_encoded_dataset_df.head(), "\n\n")
        X_hot_train, X_hot_test, y_hot_train, y_hot_test = generateTrainAndTestSplit(hot_encoded_dataset_df.values)
        
        print("""
          X: {}
          Y: {}
          TX: {}
          TY: {}
        """.format(X_hot_train.shape,  y_hot_train.shape, X_hot_test.shape, y_hot_test.shape)) 
        
        if(  row["IS_TRAINING"] == True  ):
          hot_encoded_train_features = np.append(hot_encoded_train_features, X_hot_train, axis=0)
          hot_encoded_train_labels   = np.append(hot_encoded_train_labels  , y_hot_train, axis=0)
          hot_encoded_test_features  = np.append(hot_encoded_test_features , X_hot_test , axis=0)
          hot_encoded_test_labels    = np.append(hot_encoded_test_labels   , y_hot_test , axis=0)
        else:
          print("\nAPPENDING TO VALIDATION DATA")
          hot_encoded_val_features  = np.append(hot_encoded_val_features, X_hot_train, axis=0)
          hot_encoded_val_labels    = np.append(hot_encoded_val_labels  , y_hot_train, axis=0)
          hot_encoded_val_features  = np.append(hot_encoded_val_features , X_hot_test , axis=0)
          hot_encoded_val_labels    = np.append(hot_encoded_val_labels   , y_hot_test , axis=0)

        print("\n\n", 10*"*", "FASTA TO TETRA-NUCLEOTDE FRECUENCY", 10*"*")
        tetra_n_array_positive = fastaToTetraNucletideDic(p_df.values, 1)
        tetra_n_array_negative = fastaToTetraNucletideDic(n_df.values, 0)
        joined_df              = joinPositiveAndNegative(tetra_n_array_positive , tetra_n_array_negative )
        joined_df              = joined_df.fillna(0)
        print("\nHEAD-FASTA TO TETRA-NUCLEOTDE FRECUENCY")
        print("\n\n", joined_df.head() , "\n\n" )
        X_train, X_test, y_train, y_test = generateTrainAndTestSplit(joined_df.values)
        
        print("""
          X: {}
          Y: {}
          TX: {}
          TY: {}
        """.format( X_train.shape, y_train.shape, X_test.shape, y_test.shape  ))

        if(  row["IS_TRAINING"] == True  ):
          tetra_freq_train_features = np.append(tetra_freq_train_features, X_train, axis=0)
          tetra_freq_train_labels   = np.append(tetra_freq_train_labels, y_train, axis=0)
          tetra_freq_test_features  = np.append(tetra_freq_test_features, X_test, axis=0)
          tetra_freq_test_labels    = np.append(tetra_freq_test_labels, y_test, axis=0)
        else:
          print("APPENDING TO VALIDATION DATA")
          tetra_freq_val_features = np.append(tetra_freq_val_features, X_train, axis=0)
          tetra_freq_val_labels   = np.append(tetra_freq_val_labels, y_train, axis=0)
          tetra_freq_val_features  = np.append(tetra_freq_val_features, X_test, axis=0)
          tetra_freq_val_labels    = np.append(tetra_freq_val_labels, y_test, axis=0)

        print("\n\n", 10*"*", "RNN DATA PROCESSING", 10*"*")

        tetran_word_dataset = fasta_to_tetranucleotide_list(p_df.values, n_df.values)
        tetran_word_dataset = tetran_word_dataset.dropna()
        print("\n\n", tetran_word_dataset.head(), "\n\n" )
        X_tetran_train, X_tetran_test, y_tetran_train, y_tetran_test = generateTrainAndTestSplit( tetran_word_dataset.values)
        
        print("""\n
          X:  {}
          Y:  {}
          TX: {}
          TY: {}
          COMPLETE:        {}
          COMPLETE+LABELS: {}
        """.format( 
          np.array(X_tetran_train).shape , 
          np.array(y_tetran_train).shape,
          np.array(X_tetran_test).shape , 
          np.array(y_tetran_test).shape, 
          np.array(tetran_word_dataset.iloc[:, :-1].values).shape, 
          np.array(tetran_word_dataset.values).shape 
        ))

        if(  row["IS_TRAINING"] == True  ):
          rnn_token_train_features = np.append(rnn_token_train_features, X_tetran_train, axis=0)
          rnn_token_train_labels   = np.append(rnn_token_train_labels  , y_tetran_train, axis=0)
          rnn_token_test_features  = np.append(rnn_token_test_features , X_tetran_test , axis=0)
          rnn_token_test_labels    = np.append(rnn_token_test_labels   , y_tetran_test , axis=0)
        else:
          print("APPENDING TO VALIDATION DATA")
          rnn_token_val_features = np.append(rnn_token_val_features , X_tetran_train, axis=0)
          rnn_token_val_labels   = np.append(rnn_token_val_labels   , y_tetran_train, axis=0)
          rnn_token_val_features = np.append(rnn_token_val_features , X_tetran_test , axis=0)
          rnn_token_val_labels   = np.append(rnn_token_val_labels   , y_tetran_test , axis=0)
        global_rnn_complete       = np.append(global_rnn_complete      , tetran_word_dataset.iloc[:, :-1].values , axis=0)

      except Exception as e:
        print('\n\nFAILED : \n\n'+ str(e))
        print(traceback.format_exc())
        error_msg = str(e)

    if(nextflow_path is not None):
      w.writerow([ row['BACTERIA'] ,  stdout,  stderr,  error_msg ])
    
    bar.update(i)
    bacteria_end_time = datetime.datetime.now().time().strftime('%H:%M:%S')
    bacteria_total_time=(datetime.datetime.strptime(bacteria_end_time,'%H:%M:%S') - datetime.datetime.strptime(bacteria_start_time,'%H:%M:%S'))
    print("\n\nBACTERIA: ", row['BACTERIA'] , " - TOTAL ELAPSED TIME: ",bacteria_total_time )


    
  print("\n\nTOKENIZING RNN DATASET\n\n")
  str_global_rnn_complete = tetranucleotide_list_to_string_list(global_rnn_complete)
  str_rnn_token_train_features = tetranucleotide_list_to_string_list(rnn_token_train_features)
  str_rnn_token_test_features  = tetranucleotide_list_to_string_list(rnn_token_test_features)
  str_rnn_token_val_features   = tetranucleotide_list_to_string_list(rnn_token_val_features)

  tokenizer = Tokenizer()
  tokenizer.fit_on_texts(str_global_rnn_complete)
  vocab_size = len(tokenizer.word_index)+1

  print("\nTokenizer Summary")
  print("\n document_count: ", tokenizer.document_count)
  print("\n vocab size: ", vocab_size )

  rnn_token_train_features = tokenizer.texts_to_sequences(str_rnn_token_train_features)
  rnn_token_test_features  = tokenizer.texts_to_sequences(str_rnn_token_test_features)
  rnn_token_val_features   = tokenizer.texts_to_sequences(str_rnn_token_val_features)
  # X_train_pad = pad_sequences(rnn_token_train_features, maxlen=37, padding="post")
  # X_test_pad  = pad_sequences(rnn_token_test_features, maxlen=37, padding="post")
  rnn_token_train_features = np.array(rnn_token_train_features)
  rnn_token_test_features  = np.array(rnn_token_test_features)
  rnn_token_val_features  = np.array(rnn_token_val_features)

  print(
    "\n\nTOTAL HOT ENCODING FEATURES"
    "\nHOT ENCODED FEATURE TRAIN", hot_encoded_train_features.shape,
    "\nHOT ENCODED LABELS  TRAIN", hot_encoded_train_labels.shape,
    "\nHOT ENCODED FEATURE TEST", hot_encoded_test_features.shape,
    "\nHOT ENCODED LABELS  TEST", hot_encoded_test_labels.shape,
    "\nHOT ENCODED FEATURE VAL", hot_encoded_val_features.shape,
    "\nHOT ENCODED LABELS  VAL", hot_encoded_val_labels.shape,
    "\n"
  ) 
  print(
    "\n\nTOTAL TETRA-NUCLEOTDE FRECUENCY FEATURES"
    "\nTETRA-FREQ FEATURE TRAIN", tetra_freq_train_features.shape,
    "\nTETRA-FREQ LABELS  TRAIN", tetra_freq_train_labels.shape,
    "\nTETRA-FREQ FEATURE TEST", tetra_freq_test_features.shape,
    "\nTETRA-FREQ LABELS  TEST", tetra_freq_test_labels.shape,
    "\nTETRA-FREQ FEATURE VAL", tetra_freq_val_features.shape,
    "\nTETRA-FREQ LABELS  VAL", tetra_freq_val_labels.shape,
    "\n"
  )
  print(
    "\n\nTOTAL RNN TETRANUCLEOTIDE STRING TOKEN SEQUENCES FEATURES"
    "\nRNN TOKEN FEATURE TRAIN", rnn_token_train_features.shape,
    "\nRNN TOKEN LABELS  TRAIN", rnn_token_train_labels.shape,
    "\nRNN TOKEN FEATURE TEST", rnn_token_test_features.shape,
    "\nRNN TOKEN LABELS  TEST", rnn_token_test_labels.shape,
    "\nRNN TOKEN FEATURE VAL", rnn_token_val_features.shape,
    "\nRNN TOKEN LABELS  VAL", rnn_token_val_labels.shape,
    "\nRNN TOKEN ALL", global_rnn_complete.shape,
    "\nVocab", vocab_size,
    "\n"
  )
  # Save files
  if(save_data):
      saveData(
        hot_encoded_train_features, 
        hot_encoded_train_labels, 
        hot_encoded_test_features, 
        hot_encoded_test_labels,  
        hot_encoded_val_features, 
        hot_encoded_val_labels,  
        tetra_freq_train_features, 
        tetra_freq_train_labels, 
        tetra_freq_test_features, 
        tetra_freq_test_labels, 
        tetra_freq_val_features, 
        tetra_freq_val_labels, 
        rnn_token_train_features, 
        rnn_token_train_labels,  
        rnn_token_test_features, 
        rnn_token_test_labels, 
        rnn_token_val_features, 
        rnn_token_val_labels, 
        vocab_size, tokenizer,
        save_csv
      )
      try:
        print("\n\nDeleting Temporary Files\n\n")
        os.system('rm -rf __pycache__')
        os.system('rm -rf .nextflow')
        #os.system('rm -rf outDirOriginal')
        #os.system('rm -rf work')
        #os.system('rm .nextflow.*')
        #os.system('mv -v *.genome ./data')
        #os.system('mkdir -p ./data/bacteria')
        #os.system('mv ./outDir/* ./data/bacteria')
        #os.system('rm -rf ./outDir')
      except Exception as e:
        print("\n\nError deleting temporary data. "+str(e))
  else:
    print("NOT SAVING BINARY DATA")
  
  end_time   = datetime.datetime.now().time().strftime('%H:%M:%S')
  total_time = (datetime.datetime.strptime(end_time,'%H:%M:%S') - datetime.datetime.strptime(start_time,'%H:%M:%S'))
  print("\n\nTOTAL ELAPSED TIME: ",total_time )

  return hot_encoded_train_features, \
    hot_encoded_train_labels, \
    hot_encoded_test_features, \
    hot_encoded_test_labels,  \
    hot_encoded_val_features, \
    hot_encoded_val_labels,  \
    tetra_freq_train_features, \
    tetra_freq_train_labels, \
    tetra_freq_test_features, \
    tetra_freq_test_labels, \
    tetra_freq_val_features, \
    tetra_freq_val_labels, \
    rnn_token_train_features, \
    rnn_token_train_labels,  \
    rnn_token_test_features, \
    rnn_token_test_labels, \
    rnn_token_val_features, \
    rnn_token_val_labels, \
    vocab_size, tokenizer

def loadData(data_path="./"):
  print("LOAD DATA")
  hot_encoded_train_features = pickle.load(open(data_path+'TRAIN_HOT_ENCODED_FEATURES.data', 'rb'))
  hot_encoded_train_labels   = pickle.load(open(data_path+'TRAIN_HOT_ENCODED_LABELS.data'  , 'rb'))
  hot_encoded_test_features  = pickle.load(open(data_path+'TEST_HOT_ENCODED_FEATURES.data' , 'rb'))
  hot_encoded_test_labels    = pickle.load(open(data_path+'TEST_HOT_ENCODED_LABELS.data'   , 'rb'))
  hot_encoded_val_features   = pickle.load(open(data_path+'VAL_HOT_ENCODED_FEATURES.data'  , 'rb'))
  hot_encoded_val_labels     = pickle.load(open(data_path+'VAL_HOT_ENCODED_LABELS.data'    , 'rb'))

  tetra_freq_train_features     = pickle.load(open(data_path+'TRAIN_TETRA_FREQ_FEATURES.data'    , 'rb'))
  tetra_freq_train_labels       = pickle.load(open(data_path+'TRAIN_TETRA_FREQ_LABELS.data'      , 'rb'))
  tetra_freq_test_features      = pickle.load(open(data_path+'TEST_TETRA_FREQ_FEATURES.data'     , 'rb'))
  tetra_freq_test_labels        = pickle.load(open(data_path+'TEST_TETRA_FREQ_LABELS.data'       , 'rb'))
  tetra_freq_val_features       = pickle.load( open(data_path+'VAL_TETRA_FREQ_FEATURES.data'     , 'rb'))
  tetra_freq_val_labels         = pickle.load( open(data_path+'VAL_TETRA_FREQ_LABELS.data'       , 'rb'))


  rnn_token_train_features = pickle.load(open(data_path+'TRAIN_RNN_TOKENS_FEATURES.data', 'rb'))
  rnn_token_train_labels   = pickle.load(open(data_path+'TRAIN_RNN_TOKENS_LABELS.data'  , 'rb'))
  rnn_token_test_features  = pickle.load(open(data_path+'TEST_RNN_TOKENS_FEATURES.data' , 'rb'))
  rnn_token_test_labels    = pickle.load(open(data_path+'TEST_RNN_TOKENS_LABELS.data'   , 'rb'))
  rnn_token_val_features   = pickle.load(open(data_path+'VAL_RNN_TOKENS_FEATURES.data'  , 'rb'))
  rnn_token_val_labels     = pickle.load(open(data_path+'VAL_RNN_TOKENS_LABELS.data'    , 'rb'))
  vocab_size                = pickle.load(open(data_path+'vocab_size.data'               , 'rb'))
  tokenizer                 = pickle.load(open(data_path+'tokenizer.data'                , 'rb'))

  print(
    "\n\nTOTAL HOT ENCODING FEATURES"
    "\nHOT ENCODED FEATURE TRAIN", hot_encoded_train_features.shape,
    "\nHOT ENCODED LABELS  TRAIN", hot_encoded_train_labels.shape,
    "\nHOT ENCODED FEATURE TEST", hot_encoded_test_features.shape,
    "\nHOT ENCODED LABELS  TEST", hot_encoded_test_labels.shape,
    "\nHOT ENCODED FEATURE VAL", hot_encoded_val_features.shape,
    "\nHOT ENCODED LABELS  VAL", hot_encoded_val_labels.shape,
    "\n"
  ) 
  print(
    "\n\nTOTAL TETRA-NUCLEOTDE FRECUENCY FEATURES"
    "\nTETRA-FREQ FEATURE TRAIN", tetra_freq_train_features.shape,
    "\nTETRA-FREQ LABELS  TRAIN", tetra_freq_train_labels.shape,
    "\nTETRA-FREQ FEATURE TEST", tetra_freq_test_features.shape,
    "\nTETRA-FREQ LABELS  TEST", tetra_freq_test_labels.shape,
    "\nTETRA-FREQ FEATURE VAL", tetra_freq_val_features.shape,
    "\nTETRA-FREQ LABELS  VAL", tetra_freq_val_labels.shape,
    "\n"
  )
  print(
    "\n\nTOTAL RNN TETRANUCLEOTIDE STRING TOKEN SEQUENCES FEATURES"
    "\nRNN TOKEN FEATURE TRAIN", rnn_token_train_features.shape,
    "\nRNN TOKEN LABELS  TRAIN", rnn_token_train_labels.shape,
    "\nRNN TOKEN FEATURE TEST", rnn_token_test_features.shape,
    "\nRNN TOKEN LABELS  TEST", rnn_token_test_labels.shape,
    "\nRNN TOKEN FEATURE VAL", rnn_token_val_features.shape,
    "\nVocab", vocab_size,
    "\n"
  )
  return hot_encoded_train_features, \
  hot_encoded_train_labels, \
  hot_encoded_test_features, \
  hot_encoded_test_labels,  \
  hot_encoded_val_features, \
  hot_encoded_val_labels,   \
  tetra_freq_train_features, \
  tetra_freq_train_labels, \
  tetra_freq_test_features, \
  tetra_freq_test_labels,   \
  tetra_freq_val_features,  \
  tetra_freq_val_labels,    \
  rnn_token_train_features, \
  rnn_token_train_labels,  \
  rnn_token_test_features, \
  rnn_token_test_labels,   \
  rnn_token_val_features, \
  rnn_token_val_labels,   \
  vocab_size,  tokenizer

def saveData(
  hot_encoded_train_features, 
  hot_encoded_train_labels, 
  hot_encoded_test_features, 
  hot_encoded_test_labels,  
  hot_encoded_val_features, 
  hot_encoded_val_labels,  
  tetra_freq_train_features, 
  tetra_freq_train_labels, 
  tetra_freq_test_features, 
  tetra_freq_test_labels, 
  tetra_freq_val_features, 
  tetra_freq_val_labels, 
  rnn_token_train_features, 
  rnn_token_train_labels,  
  rnn_token_test_features, 
  rnn_token_test_labels, 
  rnn_token_val_features, 
  rnn_token_val_labels, 
  vocab_size, tokenizer ,
  save_csv=False ):

  data_root = "./data/"
  if not os.path.exists(data_root):
      os.makedirs(data_root)
  print("SAVE DATA")

  pickle.dump( hot_encoded_train_features , open(data_root+'TRAIN_HOT_ENCODED_FEATURES.data', 'wb'))
  pickle.dump( hot_encoded_train_labels   , open(data_root+'TRAIN_HOT_ENCODED_LABELS.data'  , 'wb'))
  pickle.dump( hot_encoded_test_features  , open(data_root+'TEST_HOT_ENCODED_FEATURES.data' , 'wb'))
  pickle.dump( hot_encoded_test_labels    , open(data_root+'TEST_HOT_ENCODED_LABELS.data'   , 'wb'))
  pickle.dump( hot_encoded_val_features  , open(data_root+'VAL_HOT_ENCODED_FEATURES.data' , 'wb'))
  pickle.dump( hot_encoded_val_labels    , open(data_root+'VAL_HOT_ENCODED_LABELS.data'   , 'wb'))

  pickle.dump( tetra_freq_train_features     , open(data_root+'TRAIN_TETRA_FREQ_FEATURES.data'    , 'wb'))
  pickle.dump( tetra_freq_train_labels       , open(data_root+'TRAIN_TETRA_FREQ_LABELS.data'      , 'wb'))
  pickle.dump( tetra_freq_test_features      , open(data_root+'TEST_TETRA_FREQ_FEATURES.data'     , 'wb'))
  pickle.dump( tetra_freq_test_labels        , open(data_root+'TEST_TETRA_FREQ_LABELS.data'       , 'wb'))
  pickle.dump( tetra_freq_val_features      , open(data_root+'VAL_TETRA_FREQ_FEATURES.data'     , 'wb'))
  pickle.dump( tetra_freq_val_labels        , open(data_root+'VAL_TETRA_FREQ_LABELS.data'       , 'wb'))
  
  pickle.dump( rnn_token_train_features , open(data_root+'TRAIN_RNN_TOKENS_FEATURES.data', 'wb'))
  pickle.dump( rnn_token_train_labels   , open(data_root+'TRAIN_RNN_TOKENS_LABELS.data'  , 'wb'))
  pickle.dump( rnn_token_test_features  , open(data_root+'TEST_RNN_TOKENS_FEATURES.data' , 'wb'))
  pickle.dump( rnn_token_test_labels    , open(data_root+'TEST_RNN_TOKENS_LABELS.data'   , 'wb'))
  pickle.dump( rnn_token_val_features  , open(data_root+'VAL_RNN_TOKENS_FEATURES.data' , 'wb'))
  pickle.dump( rnn_token_val_labels    , open(data_root+'VAL_RNN_TOKENS_LABELS.data'   , 'wb'))

  pickle.dump( tokenizer                 , open(data_root+'tokenizer.data', 'wb'))
  pickle.dump( vocab_size                , open(data_root+'vocab_size.data'               , 'wb'))

  if(save_csv):
    pd.DataFrame(hot_encoded_train_features).to_csv(data_root+"TRAIN_HOT_ENCODED_FEATURES.csv")
    pd.DataFrame(hot_encoded_train_labels).to_csv(data_root+"TRAIN_HOT_ENCODED_LABELS.csv")
    pd.DataFrame(hot_encoded_test_features).to_csv(data_root+"TEST_HOT_ENCODED_FEATURES.csv")
    pd.DataFrame(hot_encoded_test_labels).to_csv(data_root+"TEST_HOT_ENCODED_LABELS.csv")
    pd.DataFrame(hot_encoded_val_features).to_csv(data_root+"VAL_HOT_ENCODED_FEATURES.csv")
    pd.DataFrame(hot_encoded_val_labels).to_csv(data_root+"VAL_HOT_ENCODED_LABELS.csv")

    pd.DataFrame(tetra_freq_train_features).to_csv(data_root+"TRAIN_TETRA_FREQ_FEATURES.csv")
    pd.DataFrame(tetra_freq_train_labels).to_csv(data_root+"TRAIN_TETRA_FREQ_LABELS.csv")
    pd.DataFrame(tetra_freq_test_features).to_csv(data_root+"TEST_TETRA_FREQ_FEATURES.csv")
    pd.DataFrame(tetra_freq_test_labels).to_csv(data_root+"TEST_TETRA_FREQ_LABELS.csv")
    pd.DataFrame(tetra_freq_val_features).to_csv(data_root+"VAL_TETRA_FREQ_FEATURES.csv")
    pd.DataFrame(tetra_freq_val_labels).to_csv(data_root+"VAL_TETRA_FREQ_LABELS.csv")

    pd.DataFrame(rnn_token_train_features).to_csv(data_root+"TRAIN_RNN_TOKENS_FEATURES.csv")
    pd.DataFrame(rnn_token_train_labels).to_csv(data_root+"TRAIN_RNN_TOKENS_LABELS.csv")
    pd.DataFrame(rnn_token_test_features).to_csv(data_root+"TEST_RNN_TOKENS_FEATURES.csv")
    pd.DataFrame(rnn_token_test_labels).to_csv(data_root+"TEST_RNN_TOKENS_LABELS.csv")
    pd.DataFrame(rnn_token_val_features).to_csv(data_root+"VAL_RNN_TOKENS_FEATURES.csv")
    pd.DataFrame(rnn_token_val_labels).to_csv(data_root+"VAL_RNN_TOKENS_LABELS.csv")
  
def parseData(fasta_file, bed_file, output_dir, output_name, promoter_size=40, sliding_step=1, tetranucleotide_size=4, test_with_sample=False, test_sample_size=5000, tokenizer_path=None, data_output_type=None, save_to_disk=True  ):
  start_time = time.time()
  seq_arr     = open(fasta_file, 'r').readlines() 
  seq_info    = seq_arr[0]
  seq_id      = seq_info.split(" ")[0][1:]
  all_genome  = "".join(seq_arr[1:]).replace("\n", "")
  logs_arr    = list()
  logs_arr    = h1(logs_arr,  '{} - PROMOTER PREDICTION OF GENOME {} WITH LENGTH {} USING SLIDING WINDOW FROM {}'.format(output_name, seq_id, len(seq_arr), fasta_file ) )
  
  prom_arr, inv_prom_arr, prom_tetra_arr, inv_prom_tetra_arr, prom_tetra_arr_str, inv_prom_tetra_arr_str = None, None, None, None, None, None

  prom_arr_path           = "{}/{}/{}.data".format( output_dir, output_name, "40BP_SEQUENCES" )
  inv_prom_arr_path       = "{}/{}/{}.data".format( output_dir, output_name, "40BP_SEQUENCES_INV" )
  prom_tetra_arr_path     = "{}/{}/{}.data".format( output_dir, output_name, "TETRA_NUCLEOTIDES" )
  inv_prom_tetra_arr_path = "{}/{}/{}.data".format( output_dir, output_name, "TETRA_NUCLEOTIDES_INV" )

  if os.path.exists(prom_arr_path) and os.path.exists(inv_prom_arr_path) and os.path.exists(prom_tetra_arr_path) and os.path.exists(inv_prom_tetra_arr_path):
    logs_arr = h1(logs_arr,  "(1-3/4) LOADING EXISTING SEQUENCES AND TETRANUCLEOTIDE DATA")  
    
    prom_arr           = joblib.load(prom_arr_path)
    inv_prom_arr       = joblib.load(inv_prom_arr_path)
    prom_tetra_arr     = joblib.load(prom_tetra_arr_path)
    inv_prom_tetra_arr = joblib.load(inv_prom_tetra_arr_path)

    logs_arr = h1(logs_arr, """
    SAMPLES:
    _______________________________________
    ORIGINAL : {} \n\t {} \n
    INVERSE  : {} \n\t {} \n
    -----------
    ORIGINAL JOINT : {} \n\t {} \n
    INVERSE  JOINT : {} \n\t {} \n
    """.format (  prom_arr.shape, prom_arr[-1], inv_prom_arr.shape, inv_prom_arr[-1], prom_tetra_arr.shape, prom_tetra_arr[-1], inv_prom_tetra_arr.shape, inv_prom_tetra_arr[-1]  ))
  
    logs_arr = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) )

    if(test_with_sample):
      logs_arr = h1(logs_arr, '(3.25/4) REDUCING DATASET FOR TESTING. FROM SIZE {} TO {}.'.format(prom_arr.shape, test_sample_size ))
      prom_arr           = prom_arr[:test_sample_size]
      inv_prom_arr       = inv_prom_arr[:test_sample_size]
      prom_tetra_arr     = prom_tetra_arr[:test_sample_size]
      inv_prom_tetra_arr = inv_prom_tetra_arr[:test_sample_size]

    logs_arr = h1(logs_arr,  "\t (3.5/4) CREATING JOINT TETRANUCLEOTIDE ARRAYS")  

    prom_tetra_arr_str     = tetranucleotide_list_to_string_list(prom_tetra_arr)
    inv_prom_tetra_arr_str = tetranucleotide_list_to_string_list(inv_prom_tetra_arr)

    logs_arr = h1(logs_arr, """
    SAMPLES:
    _______________________________________
    ORIGINAL : \n\t {} \n
    INVERSE  : \n\t {} \n
    """.format ( prom_tetra_arr_str[-1], inv_prom_tetra_arr_str[-1] ))
    logs_arr = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) )

  else:
    logs_arr = h1(logs_arr,  "(1/4) CUTTING GENOME SEQUENCE INTO {} BP SEQUENCES USING SLIDING WINDOW OF STEP {}. TOTAL SAMPLES: {}. SAMPLES USING SLIDING STEP: {}".format( promoter_size, sliding_step, len(all_genome)-promoter_size-1, sliding_step ))  
    prom_arr = np.array([ all_genome[i:i+promoter_size] for i in progressbar.progressbar( range( 0, len(all_genome)-promoter_size-1 , sliding_step ) ) ])   
    logs_arr = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) )

    if(test_with_sample):
      logs_arr = h1(logs_arr, '(1.5/4) REDUCING DATASET FOR TESTING. FROM SIZE {} TO {}.'.format(prom_arr.shape, test_sample_size ))
      prom_arr           = prom_arr[:test_sample_size] 
    else:
      logs_arr = h1(logs_arr, '\t DATA IS NOT FOR TESTING. RUNNING COMPLETE DATASET OF SIZE: {}. '.format( prom_arr.shape) )
    logs_arr   = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) )

    logs_arr = h1(logs_arr,  "(2/4) OBTAINING TETRANUCLEOTIDES FROM EACH 40 BP SEQUENCES"  )
    logs_arr = h1(logs_arr,  """\t RESULTING SEQUENCES ARRAY OF SIZE {}. EACH SEQUENCE WITH SIZE OF {} bp. SAMPLE:  "{}"\n""".format( len(prom_arr), len(prom_arr[-1]), prom_arr[-1] ))
    prom_tetra_arr        = np.empty( len(prom_arr) , dtype=object )
    for i_s in progressbar.progressbar( range(len(prom_arr)) ):
      sequence            = prom_arr[i_s]
      prom_tetra_arr[i_s] = np.array([ sequence[i_t:i_t+tetranucleotide_size] for i_t in range( len(sequence) - (tetranucleotide_size - 1) ) ])
    prom_tetra_arr_str    = tetranucleotide_list_to_string_list(prom_tetra_arr)
    logs_arr = h1(logs_arr, """
    SAMPLES:
    _______________________________________
    ORIGINAL: {} \n\t {} \n\n TETRANUCLEOTODES: {} \n\t {} \n
    """.format ( prom_arr.shape, prom_arr[-1], prom_tetra_arr_str.shape, prom_tetra_arr_str[-1] ))
    logs_arr = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) )

    logs_arr = h1(logs_arr, "(3/4) GENERATING INVERSE SEQUENCES " )
    inv_prom_arr              = np.array([ str(Seq(seq[::-1]).complement()) for seq in prom_arr ])
    inv_prom_tetra_arr        = np.empty( len(inv_prom_arr), dtype=object )
    for i_s in progressbar.progressbar( range(len(inv_prom_arr)) ):
      sequence                = inv_prom_arr[i_s]
      inv_prom_tetra_arr[i_s] = np.array([sequence[i_t:i_t+tetranucleotide_size] for i_t in range(len(sequence) - (tetranucleotide_size-1)  )])
    inv_prom_tetra_arr_str    = tetranucleotide_list_to_string_list(inv_prom_tetra_arr)
    logs_arr = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time) )))
    logs_arr = h1(logs_arr, """
    SAMPLES:
    _______________________________________
    ORIGINAL : {} \n\t {} \n
    INVERSE  : {} \n\t {} \n
    -----------
    ORIGINAL JOINT : {} \n\t {} \n
    INVERSE  JOINT : {} \n\t {} \n
    """.format (  prom_arr.shape, prom_arr[-1], inv_prom_arr.shape, inv_prom_arr[-1], prom_tetra_arr_str.shape, prom_tetra_arr_str[-1], inv_prom_tetra_arr_str.shape, inv_prom_tetra_arr_str[-1]  ))
    logs_arr = h1(logs_arr,  "\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time) )))

  X, X_inv = None, None 
  if(data_output_type == "RNN"):
    logs_arr = h1(logs_arr, "\t LOADING RNN TOKENIZER" )
    tokenizer               = pickle.load(open(tokenizer_path, 'rb'))
    logs_arr = h1(logs_arr, "\t TOKENIZER   : {}".format(tokenizer))

    logs_arr = h1(logs_arr, "(4/5) PARSING PROMOTER 40BP SEQUENCES TO RNN TOKEN SEQUENCES" )
    tetra_tokens     = np.array( tokenizer.texts_to_sequences(prom_tetra_arr_str) )
    inv_tetra_tokens = np.array( tokenizer.texts_to_sequences(inv_prom_tetra_arr_str) )
    X, X_inv         = tetra_tokens, inv_tetra_tokens
    logs_arr = h1(logs_arr, """
    TOKEN SAMPLES:
    _______________________________________
    ORIGINAL : {} \n\t {} \n\t {} \n
    -----------
    INVERSE  : {} \n\t {} \n\t {} \n
    """.format( X.shape, prom_arr[-1], X[-1], X_inv.shape, X_inv[-1], inv_prom_arr[-1] ))
  elif(data_output_type == "RF-HOT"):
    logs_arr = h1(logs_arr, "(4/5)  PARSING PROMOTER 40BP SEQUENCES TO HOT ENCODING FORMAT" )
    tetra_hot     = tetraToHotEncoding(prom_arr)
    inv_tetra_hot = tetraToHotEncoding(inv_prom_arr)
    X, X_inv      = tetra_hot, inv_tetra_hot
    logs_arr = h1(logs_arr, """
    HOT ENCODING SAMPLES:
    _______________________________________
    ORIGINAL : {} \n\t {} \n\t {} \n
    -----------
    INVERSE  : {} \n\t {} \n\t {} \n
    """.format( X.shape, prom_arr[-1], X[-1], X_inv.shape, X_inv[-1], inv_prom_arr[-1] ))
  elif(data_output_type == "RF-TETRA"):
    logs_arr = h1(logs_arr, "(4/5) PARSING PROMOTER-40BP SEQUENCES TO TETRA-FREQUENCIES FORMAT" )
    tetra_freq     = promoterToTetraFreq(prom_arr)
    inv_tetra_freq = promoterToTetraFreq(inv_prom_arr)
    X, X_inv = tetra_freq.values, inv_tetra_freq.values
  else:
    logs_arr = h1(logs_arr, "ERROR: NO DATA OUTPUT TYPE SPECIFIED.")
    raise "ERROR: NO DATA OUTPUT TYPE SPECIFIED."
  
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  if not os.path.exists( "{}/{}".format(output_dir, output_name) ):
    os.makedirs( "{}/{}".format(output_dir, output_name) )
  if(save_to_disk):
    logs_arr = h1(logs_arr, "(5/5) SAVING FILES" )
    
    logs_arr = h1(logs_arr, "\tSAVING X: {}".format(X.shape) )
    joblib.dump(X                 , "{}/{}/{}.data".format( output_dir, output_name, data_output_type))
    
    logs_arr = h1(logs_arr, "\tSAVING X_INV: {}".format(X_inv.shape) )
    joblib.dump(X_inv             , "{}/{}/{}_INV.data".format( output_dir, output_name, data_output_type))
  else:
    logs_arr = h1(logs_arr, "NOT SAVING TO DISK.")
  logs_arr = h1(logs_arr,  "\tTIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time) )))

  with open("{}/{}/{}_LOGS.txt".format( output_dir, output_name, data_output_type), 'w') as f:
    for item in logs_arr:
        f.write("%s\n" % item)
  return X, X_inv


def increaseDataSize(
  bacteria_name,fasta_file,  
  p_bed_path, n_bed_path, 
  genome_file_path, 
  output_folder, 
  both_sides_padding_size=20
):
  print("INCREASING BED SEQ SIZES: ", bacteria_name, "\n", 10*"__") 
  if not os.path.exists(output_folder):
    print("CREATING: ", output_folder)
    os.mkdir(output_folder)
  output_folder = "{}/{}".format(output_folder,bacteria_name)
  if not os.path.exists(output_folder):
    os.mkdir(output_folder)
    
  input_p_bed_df = BedTool(p_bed_path).to_dataframe()
  input_n_bed_df = BedTool(n_bed_path).to_dataframe()
  print("INPUT - P {} N {}".format(input_p_bed_df.shape, input_n_bed_df.shape))
  
  print("BED SLOP POSITIVE")
  os.system("bedtools slop -b {} -i {} -g {} > {}/positive.bed".format( 
    both_sides_padding_size,
    p_bed_path, 
    genome_file_path, 
    output_folder
  ))
  print("BED SLOP NEGATIVES")
  os.system("bedtools slop -b {} -i {} -g {} > {}/negative.bed".format( 
    both_sides_padding_size,
    n_bed_path, 
    genome_file_path, 
    output_folder
  ))
  print("GET FASTA POSITIVE")
  os.system("bedtools getfasta -fi {} -bed {}/positive.bed -s -fo {}/positive.fasta".format( 
    fasta_file,
    output_folder,
    output_folder 
  ))
  print("GET FASTA NEGATIVES")
  os.system("bedtools getfasta -fi {} -bed {}/negative.bed -s -fo {}/negative.fasta".format( 
    fasta_file,
    output_folder,
    output_folder 
  ))
  p_df, n_df = readData(
    "{}/positive.fasta".format(output_folder), 
    "{}/negative.fasta".format(output_folder), 
  )
  p_bed_df = BedTool("{}/positive.bed".format(output_folder)).to_dataframe()
  n_bed_df = BedTool("{}/negative.bed".format(output_folder)).to_dataframe()
  print(
    "P: ", p_df.shape, " N: ",  n_df.shape, 
    "P: ", p_bed_df.shape, " N: ", n_bed_df.shape 
  )
  p_bed_df["sequence"] = p_df.values
  n_bed_df["sequence"] = n_df.values
  p_bed_df["label"]    = [1]*len(p_df)
  n_bed_df["label"]    = [0]*len(n_df)
  dataset_df = pd.concat([p_bed_df, n_bed_df])
  print("SAVING DATASET: P {}  + N {} = {}".format( p_bed_df.shape, n_bed_df.shape, dataset_df.shape ))
  p_bed_df.to_csv(   "{}/positive.csv".format(output_folder)  )
  n_bed_df.to_csv(   "{}/negative.csv".format(output_folder)  )
  dataset_df.to_csv( "{}/dataset.csv".format( output_folder)   )
  print(dataset_df.head())

def parseSequenceList(fasta_file_path, sequence_length=40, slop_mode="middle", file_encoding='utf-8', tokenizer_path=None ):
  sequence_arr = list() 
  sequences = SeqIO.parse(open(fasta_file_path,  encoding=file_encoding),'fasta')
  for seq in sequences:
    sequence_arr.append([str(seq.seq)])
  seq_df = pd.DataFrame(sequence_arr)
  sample     = seq_df.iloc[0, 0]
  sample_len = len(sample)
  seq_df = seq_df.applymap( lambda x: str(x).upper() )
  print("INPUT FILE {} \nINPUT SHAPE {} \nSAMPLE WITH LEN {}: \n{}".format(fasta_file_path, seq_df.shape, sample_len, sample))
  if(sample_len > sequence_length):
    if( slop_mode == "beginning"):
      seq_df = seq_df.applymap( lambda x: x[0:sequence_length] )
    elif( slop_mode == "middle" ):
      cutting_len_at_both_sides = int((sample_len-sequence_length)/2)
      seq_df = seq_df.applymap( lambda x: x[cutting_len_at_both_sides:cutting_len_at_both_sides+sequence_length] )
    elif( slop_mode == "end" ):
      seq_df = seq_df.applymap( lambda x: x[-sequence_length:] )
  print("OUTPUT SAMPLE WITH LEN {}: \n{}".format( len(seq_df.iloc[0,0]), seq_df.iloc[0,0] ))
  hot_seq_df, _ = fastaToHotEncoding( seq_df, pd.DataFrame() )
  tetra_freq_df = fastaToTetraNucletideDic( seq_df.values, 1 )
  
  rnn_token_df = None
  if( tokenizer_path != None ):
    tokenizer          = pickle.load(open( tokenizer_path, 'rb' ))
    tetra_seq_list     = fasta_to_tetranucleotide_list(seq_df.values, np.array([]) )
    tetra_seq_list_str = tetranucleotide_list_to_string_list( tetra_seq_list.iloc[:, :-1].values )
    rnn_token_df = tokenizer.texts_to_sequences(tetra_seq_list_str)
  return seq_df, hot_seq_df, tetra_freq_df, pd.DataFrame(rnn_token_df)

  
 
# if __name__ == "__main__":
#   val_bacteria = pd.DataFrame([ #VALIDATION
# 	  ["MYCOBACTER"   , "NC_008596.1" , "../data/bacteria/MYCOBACTER/promoter_bed_file.bed"    , "./bacteria/MYCOBACTER/MYCOBACTER.fasta"       ],
#     ["CLOSTRIDIUM"  , "NC_010001.1" , "../data/bacteria/CLOSTRIDIUM/promoter_bed_file.bed"   , "./bacteria/CLOSTRIDIUM/CLOSTRIDIUM.fasta"     ],
#     ["RHODOBACTER_1", "NC_014034.1" , "../data/bacteria/RHODOBACTER_1/promoter_bed_file.bed" , "./bacteria/RHODOBACTER_1/RHODOBACTER_1.fasta" ],
#     ["RHODOBACTER_2", "NC_014034.1" , "../data/bacteria/RHODOBACTER_2/promoter_bed_file.bed" , "./bacteria/RHODOBACTER_2/RHODOBACTER_2.fasta" ],
#     ["ECOLI"        , "NC_000913.3" , "../data/bacteria/ECOLI/promoter_bed_file.bed"         , "./bacteria/ECOLI/ECOLI.fasta"                 ],
#     ["BACILLUS"     , "CP002927.1"  , "../data/bacteria/BACILLUS/promoter_bed_file.bed"      , "./bacteria/BACILLUS/BACILLUS.fasta"           ],
#     # ["AERUGINOSA"   , "NC_002516.2" , "../data/bacteria/AERUGINOSA/promoter_bed_file.bed"    , "../bacteria/AERUGINOSA/AERUGINOSA.fasta"       ],
  
#     ["SPYOGENE"   , "LR031521.1" , "../data/bacteria/SPYOGENE/promoter_bed_file.bed"    , "./bacteria/SPYOGENE/SPYOGENE.fasta"       ],
#     ["STYPHIRMURIUM"   , "NC_016810.1" , "../data/bacteria/STYPHIRMURIUM/promoter_bed_file.bed"    , "./bacteria/STYPHIRMURIUM/STYPHIRMURIUM.fasta"       ],
#     ["CPNEUMONIAE"   , "NC_000922.1" , "../data/bacteria/CPNEUMONIAE/promoter_bed_file.bed"    , "./bacteria/CPNEUMONIAE/CPNEUMONIAE.fasta"       ],
#     ["SONEIDENSIS"   , "NC_004347.2" , "../data/bacteria/SONEIDENSIS/promoter_bed_file.bed"    , "./bacteria/SONEIDENSIS/SONEIDENSIS.fasta"       ],
#     ["LINTERROGANS"   , "NZ_LT962963.1" , "../data/bacteria/LINTERROGANS/promoter_bed_file.bed"    , "./bacteria/LINTERROGANS/LINTERROGANS.fasta"       ],
#     ["SCOELICOLOR"   , "NC_003888.3" , "../data/bacteria/SCOELICOLOR/promoter_bed_file.bed"    , "./bacteria/SCOELICOLOR/SCOELICOLOR.fasta"       ],
#  	], columns=["NAME", "ID", "BED_PATH", "FASTA_PATH" ])

#   for i_b, bacteria in val_bacteria.iterrows() :
#     increaseDataSize(
#       bacteria_name=bacteria["NAME"],
#       fasta_file=bacteria["FASTA_PATH"], 
#       p_bed_path="./data/bacteria_1_1/{}/positive.bed".format( bacteria["NAME"] ),
#       n_bed_path="./data/bacteria_1_1/{}/negative.bed".format( bacteria["NAME"] ),
#       genome_file_path="./data/{}.genome".format( bacteria["NAME"] ),
#       output_folder="250_bacteria_1_1",
#       both_sides_padding_size=106 #531# FOR G4PROMFINDER & BPROM 20 = 20 + 40 + 20 = 80, BTSSFINDER 106 = 106 + 40 +106 = 252, 531+40+531 = 1102
#     )

# if __name__ == "__main__":
#   generateData(
#     bacteria_index= None,
#     save_csv=True, 
#     save_data=True, 
#     outDir="./outDir/", 
#     nextflow_path= './nextflow',
#     nextflow_pipeline = "pipeline_unbalanced_without_docker.nf", # 'pipeline_without_docker.nf' 'pipeline_without_docker.nf',# 
#     manually_balance_data = False
#   )
