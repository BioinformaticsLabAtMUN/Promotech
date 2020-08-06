import argparse
import sys, os
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtWidgets, uic
from ui.GUI import Promotech_UI
from genome.process_genome import genome40NTSequencesToHotEncoding, predictGenomeSequences
# export DISPLAY=:0.0

model_type = "RF-HOT"


def load_model(type):
  pass


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-V", "--version"         , action='version'    , version='PROMOTECH V1.0')
  parser.add_argument("-GUI", "--gui"           , help="Show interactive GUI."                 , action="store_true" )
  parser.add_argument("-P", "--predict-sequence", help="Predict 40 nucleotides FASTA sequence.", action="store_true" )
  parser.add_argument("-PG", "--parse-genome"   , help="Parse Whole Genome Before Prediction." , action="store_true" )
  parser.add_argument("-G",  "--predict-genome" , help="Predict entire genome in a FASTA sequence. Make sure to have used" , action="store_true" )
  parser.add_argument("-F", "--fasta"    , help="FASTA sequences file. ", nargs=1, default=None) #, type=argparse.FileType('r')
  parser.add_argument("-M", "--model"          , help='Type of model used. ["RF-HOT", "RF-TETRA", "GRU", "LSTM"]', choices=["RF-HOT", "RF-TETRA", "GRU", "LSTM"], default="RF-HOT")
  parser.add_argument("-TS", "--test-samples"  , help="", default=None, type=int)
  # parser.add_argument( "-RT", "--retrain"  , help="Retrain a model. " , action="store_true"  )

  args = parser.parse_args()
  if args.gui:
    app = QApplication([])
    widget = Promotech_UI(ui_path="ui/form.ui")
    widget.show()
    sys.exit(app.exec_())
  else:
    print("""
          PROMOTECH
          MODE         : {}
          ML MODEL     : {}
          INPUT TYPE   : {}
          INPUT        : {}
          TEST SAMPLES : {}
          """.format(
            "INTERACTIVE GUI" if args.gui else "COMMAND-LINE",
            "FASTA FILE" if args.fasta else "40NT SEQUENCE",
            args.model,
            args.fasta,
            args.test_samples,
    ))
    if args.predict:
      pass
    elif(args.parse_genome):
      fasta_file_path = args.fasta
      if fasta_file_path is None:
        raise ValueError("Argument (--fasta, -F) is missing.")
      # clear && python promotech.py -G -PG -TS 20000 -F tests/genome/ECOLI_2.fasta
      genome40NTSequencesToHotEncoding(
        fasta_file_path=args.fasta[0], 
        out_dir="results",
        test_sample_size=args.test_samples,
      )
    elif args.predict_genome:
      # clear && python promotech.py -G 
      predictGenomeSequences(
        out_dir    ="results",
        model_type = args.model,
        model_file_path= os.path.join("models", "{}.model".format(args.model) ) ,
      )


# USAR PATH JOIN PARA UNIR OUT_DIR CON DEMAS
