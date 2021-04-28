
import pandas as pd, numpy as np, time, progressbar
from pathlib import Path
import joblib, traceback, os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
from Bio import SeqIO

sys.path.append("{}/../core".format(dir_path))

from core.database import Database


class iPromoter2LDatabase(Database):
    def __init__(self, filepath, extra_args={}):
        Database.__init__(self, filepath)
        self.extra_args = extra_args

    def build(self):
        if os.path.exists(self.filepath):
            self.parsed_data = SeqIO.parse(open(self.filepath), 'fasta')
            df = pd.DataFrame(columns=['id', 'seq', "size"] + list(self.extra_args.keys()))
            for seq_record in self.parsed_data: 
                dict_data = dict({"id": seq_record.id, "seq": str(seq_record.seq), "size": len(str(seq_record.seq))})
                dict_data.update(self.extra_args)
                df = df.append(dict_data, ignore_index=True)
            return df
        else:
            raise Exception("FILE: {} NOT EXISTS. DATAFRAME FROM FASTA CANNOT BE CREATED. ".format(self.filepath))
