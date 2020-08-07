import pandas as pd, numpy as np, time, progressbar
from io import StringIO
from Bio import SeqIO, Seq 
from pathlib import Path
import joblib, traceback, os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
# print("ADDING PATH: ", "{}/../core".format(dir_path) )
sys.path.append("{}/../core".format(dir_path))
from utils import fastaToHotEncodingSequences, charToBinary, print_fn

# pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 100)

def genomeSlidingWindow(fasta_file_path, log_file=None, promoter_size=40, step_size=1, print_fn=print):
  if log_file is None:
    print("NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.")
    print_fn = print
    
  if(fasta_file_path is None):
    raise "VARIABLE fasta_file_path is None. Please specify the genome FASTA file.".format()
  start_time        = time.time()
  seqs_file         = open(fasta_file_path, "r")
  seqs_file_content = seqs_file.read()
  seqs_file_str     = StringIO(seqs_file_content)
  parsed_seqs       = SeqIO.parse(seqs_file_str, "fasta")
  seqs_data         = np.array([ {"id": s.id, "seq": str(s.seq)} for s in parsed_seqs ])
  chroms            = [ s["id"]  for s in seqs_data ]
  genomes           = [ s["seq"] for s in seqs_data ]
  # CLEANING UP TO SAVE MEMORY
  if seqs_file_str is not None:
    del seqs_file_str
  if seqs_data is not None:
    del seqs_data
    
  print_fn("\n\n PRINTING CONTENT".format(), log_file) 
  for i_s, s in enumerate(genomes):
    print_fn("{}. GENOME: {} - LENGTH: {}".format(i_s+1, chroms[i_s], len(s) ), log_file)
  
  print_fn("\n\n JOINING ALL CHROMS AND SEQS INTO A SINGLE FOR TETRA-NUCLEOTIDE SLIDING WINDOW".format(), log_file) 
  chrom    = ",".join(chroms)
  genome   = "".join(genomes)
  # CLEANING UP TO SAVE MEMORY
  if chroms is not None:
    del chroms
  if genomes is not None:
    del genomes
  print_fn("\n\n JOINED GENOME: {} - LENGTH: {:,}".format( chrom, len(genome) ), log_file) 
  
  n_samples   = len( genome ) - promoter_size - 1 
  print_fn("\n\n GENERATING PROMOTER SEQUENCES WITH WINDOW-SIZE: {} AND STEP: {}. EXPECTED SAMPLES: {:,}\n\n".format( promoter_size, step_size, n_samples ), log_file) 
  cutted_seqs = np.array([ genome[i:i+promoter_size] for i in progressbar.progressbar( range( 0, n_samples, step_size ) ) ]) 
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  # CLEANING UP TO SAVE MEMORY
  if genome is not None:
    del genome
  print_fn('\n\n CUTTED {} NT SEQUENCES GENERATED SUCCESSFULLY. # OF SAMPLES: {:,} = {}. \n\tSAMPLE #1: {} \n\tSAMPLE #2: {} '.format(                                                                                                                  
    promoter_size,
    cutted_seqs.shape[0], cutted_seqs.shape, 
    cutted_seqs[0] if len(cutted_seqs)>0 else None, 
    cutted_seqs[1] if len(cutted_seqs)>1 else None 
  ), log_file)
  return chrom, cutted_seqs

def genome40NTSequencesToHotEncoding(fasta_file_path, out_dir="results", promoter_size=40, step_size=1, test_sample_size=None, print_fn=print_fn):
  if print_fn is None:
    print("NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.")
    print_fn = print
  
  if(out_dir is None):
    raise "VARIABLE out_dir is None. Please specify the output folder.".format()
  if not os.path.exists(fasta_file_path)  :
    raise ValueError("FASTA PATH {} DOES NOT EXISTS.".format(fasta_file_path))

  parent_folder      = Path(fasta_file_path).stem
  log_file           = os.path.join(out_dir, "parse.log.txt")
  start_time         = time.time()
  os.makedirs(  out_dir , exist_ok=True )
  print_fn("\n\n CREATING OUTPUT FOLDER: {}".format( out_dir ), log_file)
  
  chrom, cutted_seqs = genomeSlidingWindow(fasta_file_path=fasta_file_path, log_file=log_file, promoter_size=promoter_size, step_size=step_size, print_fn=print_fn)
  
  if test_sample_size is not None:
    print_fn("\n\n TEST SAMPLE OF SIZE: {:,} REQUESTED. REDUCING SEQUENCES FROM {:,} TO {:,}".format(  test_sample_size , cutted_seqs.shape[0] , test_sample_size ), log_file) 
    cutted_seqs = cutted_seqs[:test_sample_size]
  
  print_fn("\n\n CONVERTING {} CUTTED {} NT SEQUENCES TO HOT ENCODED SEQUENCES USING MAPPING VALUES \n\n\t{}".format( len(cutted_seqs), promoter_size, [ { nt : charToBinary(nt) } for nt in "AGCT" ] ), log_file) 
  hot_enc_seqs = fastaToHotEncodingSequences( cutted_seqs )
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  print_fn("\n\n HOT ENCODED SEQUENCES GENERATED SUCCESSFULLY. OUTPUT DATAFRAME SHAPE: {}".format( hot_enc_seqs.shape ), log_file)
  print_fn("\n\n SAMPLE: \n\n{}".format( hot_enc_seqs.head() ), log_file)
  
  forward_strand_output_path = os.path.join( out_dir, "RF-HOT.data" ) 
  chrom_output_path          = os.path.join( out_dir, "CHROM.data" )  
  seqs_output_path           = os.path.join( out_dir, "SEQS.data" )  
  print_fn("\n\n SAVING FORWARD STRAND HOT-ENCODED SEQUENCES TO BINARY FILE USING JOBLIB TO: {} ".format( forward_strand_output_path ), log_file)
  joblib.dump(  chrom        ,  chrom_output_path )
  joblib.dump(  cutted_seqs  , seqs_output_path   )
  joblib.dump(  hot_enc_seqs , forward_strand_output_path )
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  print_fn("\n\n FILE SAVED SUCCESSFULLY AT: \n\t{}".format( forward_strand_output_path ), log_file)
  
  print_fn("\n\n GENERATING INVERSE STRAND SEQUENCES. ".format( ), log_file)
  inv_cutted_seqs    = np.array([ str(Seq.Seq(cutted_seqs[i_s][::-1]).complement()) for i_s in progressbar.progressbar( range( 0, len(cutted_seqs) ) )  ])
  print_fn('\n\n INVERSE STRAND SEQUENCES GENERATED SUCCESSFULLY. # OF SAMPLES: {:,}. \n\tSAMPLE: \n\t\tORIGINAL : {} \n\t\tINVERSE  : {}'.format(                                                                                                                  
    inv_cutted_seqs.shape[0], 
    cutted_seqs[0]     , 
    inv_cutted_seqs[0] ,
  ), log_file)
  inv_seqs_output_path           = os.path.join( out_dir, "SEQS-INV.data" )  
  joblib.dump(  inv_cutted_seqs  , inv_seqs_output_path   )
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  
  print_fn("\n\n CONVERTING {} INVERSE STRAND {} NT SEQUENCES TO HOT ENCODED SEQUENCES USING MAPPING VALUES \n\n\t{}".format( len(inv_cutted_seqs), promoter_size, [ { nt : charToBinary(nt) } for nt in "AGCT" ] ), log_file) 
  hot_enc_inv_seqs = fastaToHotEncodingSequences( inv_cutted_seqs )
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ), log_file )
  print_fn("\n\n HOT ENCODED SEQUENCES GENERATED SUCCESSFULLY. OUTPUT DATAFRAME SHAPE: {}".format( hot_enc_inv_seqs.shape ), log_file)
  print_fn("\n\n SAMPLE: \n\n{}".format( hot_enc_inv_seqs.head() ), log_file)
  
  inverse_strand_output_path = os.path.join( out_dir, "RF-HOT-INV.data" ) 
  print_fn("\n\n SAVING INVERSE STRAND SEQUENCES TO BINARY FILE USING JOBLIB TO: {} ".format( inverse_strand_output_path ), log_file)
  joblib.dump( hot_enc_inv_seqs, inverse_strand_output_path )
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)
  print_fn("\n\n FILE SAVED SUCCESSFULLY AT: \n\t{}".format( inverse_strand_output_path ), log_file)
  
  return hot_enc_seqs
    

def predictGenomeSequences(
  model_file_path,
  out_dir="results",
  model_type="RF-HOT",
  threshold=0.5,
  print_fn=print_fn,
):
  if print_fn is None:
    print("NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.")
    print_fn = print
    
  log_file              = os.path.join(out_dir, "predict.log.txt")
  chrom_output_path     = os.path.join( out_dir, "CHROM.data" )  
  seqs_output_path      = os.path.join( out_dir, "SEQS.data" )  
  inv_seqs_output_path  = os.path.join( out_dir, "SEQS-INV.data" )  
  forward_strand_hot_enc_seqs_file= os.path.join("results", "{}.data".format( model_type ) ) 
  inverse_strand_hot_enc_seqs_file= os.path.join("results", "{}-INV.data".format( model_type ) ) 
  
  if not os.path.exists(out_dir)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE PARSE THE GENOME FILE FIRST.".format(out_dir))
  if not os.path.exists(forward_strand_hot_enc_seqs_file)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE PARSE THE GENOME FILE FIRST.".format(forward_strand_hot_enc_seqs_file))
  if not os.path.exists(inverse_strand_hot_enc_seqs_file)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE PARSE THE GENOME FILE FIRST.".format(inverse_strand_hot_enc_seqs_file))
  if not os.path.exists(chrom_output_path)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE PARSE THE GENOME FILE FIRST.".format(chrom_output_path))
  if not os.path.exists(seqs_output_path)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE PARSE THE GENOME FILE FIRST.".format(seqs_output_path))
  if not os.path.exists(inv_seqs_output_path)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE PARSE THE GENOME FILE FIRST.".format(inv_seqs_output_path))
  if not os.path.exists(model_file_path)  :
    raise ValueError("FILE PATH {} DOES NOT EXISTS. PLEASE MAKE SURE TO ADD YOUR MODEL '{}'(.model or .h5) FILE TO THE 'models' FOLDER.".format(model_file_path, model_type))
  
  start_time        = time.time()
  print_fn("\n\n LOADING FORWARD STRAND SEQUENCES CONVERTED TO HOT-ENCODED SEQUENCES: {} WITH SIZE: {:,.2f} MB".format( 
    forward_strand_hot_enc_seqs_file, 
    Path(forward_strand_hot_enc_seqs_file).stat().st_size  / 1000000
  ), log_file) 
  X     = joblib.load(forward_strand_hot_enc_seqs_file)
  print_fn("\n\n LOADING FORWARD STRAND SEQUENCES CONVERTED TO HOT-ENCODED SEQUENCES AT: {} WITH SIZE: {:,.2f} MB".format( 
    inverse_strand_hot_enc_seqs_file ,
    Path(inverse_strand_hot_enc_seqs_file).stat().st_size  / 1000000
  ), log_file) 
  X_INV = joblib.load(inverse_strand_hot_enc_seqs_file)
  print_fn("\n\n LOADING MACHINE LEARNING MODEL AT: {} WITH SIZE: {:,.2f} MB".format( 
    model_file_path ,
    Path(model_file_path).stat().st_size  / 1000000
  ), log_file) 
  model = joblib.load(model_file_path)
  print_fn("\n\n FORWARD STRAND SEQS: {} \nINVERSE STRAND SEQS: {} \nML-MODEL: \n\n{}".format(X.shape, X_INV.shape, str(model)), log_file)
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)  

  print_fn("\n\n GENERATING PREDICTIONS FOR FORWARD STRAND SEQUENCES WITH SHAPE: {}".format(X.shape), log_file)
  y_probs = model.predict_proba(X)
  y_pred  = y_probs[:, 1] if y_probs.shape[1] == 2 else y_probs[:, 0]
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)  
  
  print_fn("\n\n GENERATING PREDICTIONS FOR INVERSE STRAND SEQUENCES WITH SHAPE: {}".format(X.shape), log_file)
  y_inv_probs = model.predict_proba(X_INV)
  y_inv_pred  = y_inv_probs[:, 1] if y_inv_probs.shape[1] == 2 else y_inv_probs[:, 0]
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)  
  
  print_fn("\n\n PREDICTIONS GENERATED SUCCESSFULLY".format( ), log_file)
  print_fn("\n\n FORWARD SEQUENCE SAMPLE: \n\tSEQ: \n\n{} \n\tPREDICTION: {:.4f}\n".format(
    "".join(map(str, X.iloc[0].values)) ,
    y_pred[0] 
  ), log_file)
  print_fn("\n\n INVERSE SEQUENCE SAMPLE: \n\tSEQ: {} \n\tPREDICTION: {:.4f}\n".format(
    "".join(map(str, X_INV.iloc[0].values)),
    y_inv_pred[0] 
  ), log_file)
  
  print_fn("\n\n FORWARD STRAND PREDICTIONS ABOVE THRESHOLD: {:,} and BELOW: {:,} FROM TOTAL {:,}".format( len(y_pred[y_pred >= threshold]), len(y_pred[y_pred < threshold]), len(y_pred) ), log_file)
  print_fn("\n\n INVERSE STRAND PREDICTIONS ABOVE THRESHOLD: {:,} and BELOW: {:,} FROM TOTAL {:,}".format( len(y_inv_pred[y_inv_pred >= threshold]), len(y_inv_pred[y_inv_pred < threshold]), len(y_inv_pred) ), log_file)
  
  chrom                 = joblib.load( chrom_output_path )
  seqs                  = joblib.load( seqs_output_path )
  inv_seqs              = joblib.load( inv_seqs_output_path )
  print_fn("\n\n GENERATING DETECTED PROMOTERS' BED FILE BASED ON THRESHOLD: {:.3f} FOR CHROM: {}. # SEQS: {} # INV SEQS: {}. # SEQS ABOVE THRES: {} # INV SEQS ABOVE THRES: {}".format( threshold , chrom, len(seqs), len(inv_seqs), len(y_pred[y_pred >= threshold]), len(y_inv_pred[y_inv_pred >= threshold]) ), log_file)
  
  df = pd.DataFrame(columns=['chrom', 'start', 'end', 'score', 'strand', 'sequence'])
  for i_s, s in enumerate( seqs ):
    pred_score = y_pred[i_s]
    if pred_score > threshold:
      df = df.append({ 
        'chrom': chrom, 'start': i_s, 'end': i_s+39, 
        'score': np.round(pred_score, 5), 'strand': "+", 
        'sequence': seqs[i_s] 
      }, ignore_index=True)
    inv_pred_score = y_inv_pred[i_s]
    if inv_pred_score > threshold:
      df = df.append({ 
        'chrom': chrom, 'start': i_s, 'end': i_s+39, 
        'score': np.round(inv_pred_score, 5), 'strand': "-", 
        'sequence': inv_seqs[i_s] 
      }, ignore_index=True)
  pred_file_path = os.path.join(out_dir, "genome_predictions.csv")
  print_fn("\n\n SAVING BED FILE WITH SHAPE {} TO {}. SAMPLE: \n\n{}".format( df.shape, pred_file_path, df.head() ), log_file)
  df.to_csv(pred_file_path, index=None, sep='\t', columns=None)
  print_fn("\n\t TIME ELAPSED FROM START (HOUR:MIN:SEC): {}".format( time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)) ) , log_file)  
  
  
    
