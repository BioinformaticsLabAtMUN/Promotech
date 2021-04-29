import argparse
import sys, os, numpy as np, pandas as pd
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtWidgets, uic
from ui.GUI import Promotech_UI
from genome.process_genome import parseGenome40NTSequences, predictGenomeSequences
from sequences.process_sequences import predictSequences
from benchmark.process_benchmark import run_benchmark

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-v"  , "--version"          , action='version'    , version='PROMOTECH V1.0')
  parser.add_argument("-G"  , "--gui"              , help="Show interactive GUI."                 , action="store_true" )
  parser.add_argument("-s"  , "--predict-sequences", help="Predict 40 nucleotides FASTA sequence.", action="store_true" )
  parser.add_argument("-pg" , "--parse-genome"     , help="Parse Whole Genome Before Prediction." , action="store_true" )
  parser.add_argument("-ts" , "--test-samples"     , help="Parse a limited number of sequences. This argument is used together with the -PG, --parse-genome argument.", default=None, type=int)
  parser.add_argument("-g"  ,  "--predict-genome"  , help="Predict entire genome in a FASTA sequence. Make sure to have used" , action="store_true" )
  parser.add_argument("-f"  , "--fasta"            , help="FASTA sequences file. ", nargs=1, default=None) #, type=argparse.FileType('r')
  parser.add_argument("-m"  , "--model"            , help='Type of model used. The options are: ["RF-HOT", "RF-TETRA", "GRU", "LSTM"]', choices=["RF-HOT", "RF-TETRA", "GRU", "LSTM"], default="RF-HOT")
  parser.add_argument("-t"  , "--threshold"        , help='Prediction threshold.', type=float, default=0.5)
  parser.add_argument("-o"  , "--output-dir"       , help='Output folder.', type=str, default="RF-HOT")
  parser.add_argument("-i"  , "--input-dir"        , help='Input folder.', type=str, default="RF-HOT")

  parser.add_argument("-b"  , "--benchmark"        , help='Run Benchmark', action="store_true")
  parser.add_argument("-bm" , "--benchmark-model"  , help='Model to use during the benchmark. The options are: ["iPromoter2L"]', choices=["iPromoter2L"], default="iPromoter2L")
  # parser.add_argument( "-RT", "--retrain"  , help="Retrain a model. " , action="store_true"  )

  args = parser.parse_args()
  if args.gui:
    app = QApplication([])
    widget = Promotech_UI(
      ui_path="ui/form.ui"
      # init_function=None,
      # preprocess_seqs_fn=demo_preprocess,
      # predict_seqs_fn=demo_predict_seqs,
      # predict_gen_fn=None
      )
    widget.show()
    sys.exit(app.exec_())

  elif args.benchmark:
      run_benchmark(
        model_tag   = args.benchmark_model,
        chlg_model  = args.model,
        output_dir  = args.output_dir,
      )

  else:
    fasta_file_path = args.fasta
    print("""
          PROMOTECH
          MODE         : {}
          ML MODEL     : {}
          INPUT TYPE   : {}
          INPUT        : {}
          TEST SAMPLES : {}
          """.format(
            "INTERACTIVE GUI" if args.gui else "COMMAND-LINE",
            "FASTA FILE" if fasta_file_path else "40NT SEQUENCE",
            args.model,
            fasta_file_path,
            args.test_samples,
    ))
    if args.predict_sequences:
      if fasta_file_path is None:
        raise ValueError("Argument (--fasta, -F) is missing.")
      #  clear && python promotech.py -s -f examples/sequences/test.fasta -m "RF-HOT" -t 0.6
      predictSequences(
        fasta_file_path = fasta_file_path[0],
        out_dir         = args.output_dir,
        threshold       = args.threshold, 
        model_type      = args.model
      )
 
    elif(args.parse_genome):
      if fasta_file_path is None:
        raise ValueError("Argument (--fasta, -F) is missing.")
      # clear && python promotech.py-pg -ts 50000 -f examples/genome/ECOLI_2.fasta
      parseGenome40NTSequences(
        fasta_file_path  = fasta_file_path[0], 
        out_dir          = args.output_dir,
        test_sample_size = args.test_samples,
        data_type        = args.model,
      )
    elif args.predict_genome:
      # clear && python promotech.py -g -t 0.6 -m RF-HOT
      predictGenomeSequences(
        input_dir       = args.input_dir,
        out_dir         = args.output_dir,
        model_type      = args.model,
        threshold       = args.threshold
      )
      
