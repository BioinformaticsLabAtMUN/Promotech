import sys, os, numpy as np, pandas as pd, time, joblib, pickle, tensorflow as tf
from io import StringIO
from Bio import SeqIO
import progressbar
import traceback
dir_path = os.path.dirname(os.path.realpath(__file__))
# print("ADDING PATH: ", "{}/../core".format(dir_path) )
sys.path.append("{}/../core".format(dir_path))
from utils import fastaToHotEncodingSequences, charToBinary, print_fn, fastaToTetranucleotides, tetranucleotideToStringSentences, fastaToTetranucletideTable, dataConverter

def parseSequences(seqs, sequence_length=40, tokenizer_path="./models/tokenizer.data", data_type="RF-HOT", log_file=None, print_fn=print_fn ):
    if( len(seqs) == 0 ):
      raise ValueError("NO SEQUENCES TO PARSE. " + str(seqs) )
      return
    if log_file is None:
      print_fn = print
      print_fn("\n\n NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.", log_file)
    sample, sample_len     = seqs[0], len(seqs[0])
    seqs = np.array([ s.upper() for s in seqs ])
    print_fn("\n\n INPUT SHAPE {} \nSAMPLE WITH LEN {}: \n{}".format( len(seqs), sample_len, sample), log_file)
    if(sample_len > sequence_length):
        seqs = np.array([  s[0:sequence_length] for s in seqs ])
        print_fn("\n\n OUTPUT SAMPLE WITH LEN {}: \n{}".format( seqs.shape, seqs[0] ), log_file)
    print_fn("\n\n CONVERTING DATA".format(), log_file)
    data_df = dataConverter( seqs=seqs, data_type=data_type, tokenizer_path=tokenizer_path, print_fn=print_fn, log_file=log_file )
    return data_df

def predictSequencesFromString(sequences_str, print_fn=print_fn, out_dir="RF-HOT", threshold=0.5, model_type="RF-HOT" , log_file=None, start_time=time.time()):
  if log_file is None:
    print_fn = print
    print_fn("\n\n NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.", log_file)
  seqs_str          = StringIO( sequences_str )
  parsed_seqs       = SeqIO.parse( seqs_str, "fasta" )
  seqs_data         = np.array([ {"id": s.id, "seq": str(s.seq)} for s in parsed_seqs ])
  print_fn("\n\n SAMPLE: \n\n{}".format( seqs_data[0] if len(seqs_data) > 0 else seqs_data ), log_file) 
  chroms            = [ s["id"]  for s in seqs_data ]
  seqs              = [ s["seq"] for s in seqs_data ]
  print_fn("\n\n # SEQS: {}. SAMPLE: {}".format( len(seqs) , seqs[0] ), log_file) 
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  
  print_fn("\n\n CONVERTING SEQUENCES TO {} DATA TYPE".format(model_type), log_file) 
  X = parseSequences(seqs, print_fn=print_fn, data_type=model_type, log_file=log_file)
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  
  model = None
  if model_type == "RF-HOT" or model_type == "RF-TETRA":
    model_path = os.path.join("models", "{}.model".format(model_type))
    print_fn("\n\n LOADING ML MODEL {}".format(model_path), log_file) 
    model = joblib.load(model_path)
  if model_type == "GRU" or model_type == "LSTM":
    model_version = "0" if model_type == "GRU" else "3"
    model_path = os.path.join("models", "{}-{}.h5".format(model_type, model_version ))
    print_fn("\n\n LOADING ML MODEL {}".format(model_path), log_file) 
    model = tf.keras.models.load_model(model_path)
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  
  print_fn("\n\n PREDICTING SEQUNCES USING: \n\n{}".format(model), log_file) 
  y_probs = model.predict_proba( X )
  y_pred  = y_probs[:, 1] if y_probs.shape[1] == 2 else y_probs[:, 0]
  df      = pd.DataFrame({"CHROM": chroms, "SEQ": seqs, "PRED": y_pred})
  pred_file_path = os.path.join(out_dir, "sequences_predictions.csv")
  print_fn("\n\t PREDICTIONS GENERATED SUCCESSFULLY. SAMPLE: \n\n{}. \n\nSAVED AT {}".format( df.head(), pred_file_path ) , log_file)
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  df.to_csv(pred_file_path, index=None, sep='\t', columns=None)
  return chroms, seqs, y_pred, X
  
def predictSequences(fasta_file_path, print_fn=print_fn, out_dir="RF-HOT", threshold=0.5, model_type="RF-HOT" ):
  start_time        = time.time()
  log_file           = os.path.join(out_dir, "sequences.log.txt")
  if print_fn is None:
    print("NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.")
    print_fn = print
  if not os.path.exists(fasta_file_path)  :
    raise ValueError("FASTA PATH {} DOES NOT EXISTS.".format(fasta_file_path))
    
  os.makedirs( out_dir , exist_ok=True )
  
  print_fn("\n\n READING FASTA FILE: {}".format(fasta_file_path), log_file) 
  seqs_file         = open( fasta_file_path, "r" )
  seqs_file_content = seqs_file.read()
  
  predictSequencesFromString(
    sequences_str=seqs_file_content, 
    print_fn=print_fn, 
    out_dir=out_dir, 
    threshold=threshold, 
    model_type=model_type ,
    log_file=log_file,
    start_time=start_time
  )
  
  
  
  
  
  