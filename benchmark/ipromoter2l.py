from benchmark.ipromoter2ldatabase import iPromoter2LDatabase
import os
from os import path
import pandas as pd
import numpy as np
import joblib 
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("{}/../core".format(dir_path))
from utils import dataConverter
import progressbar
import logging
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_auc_score, roc_curve
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_auc_score, roc_curve
from sklearn.metrics import matthews_corrcoef, accuracy_score, f1_score
import json

os.environ['TF_CPP_MIN_VLOG_LEVEL'] = '3'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
logging.getLogger('tensorflow').disabled = True

import warnings
warnings.filterwarnings("ignore")

class iPromoter2L:
    NON_PROMOTERS_TAG = "nonPromoter"
    def __init__(self, chlg_model, output_dir="RF-HOT"):
        self.root_path = "datasets/benchmark/iPromoter-2L"
        self.databases = {
            "sigma24": iPromoter2LDatabase(path.join(self.root_path    , "sigma24.fasta")    , {"tag": "sigma24"    , "label": 1}),
            "sigma28": iPromoter2LDatabase(path.join(self.root_path    , "sigma28.fasta")    , {"tag": "sigma28"    , "label": 1}),
            "sigma32": iPromoter2LDatabase(path.join(self.root_path    , "sigma32.fasta")    , {"tag": "sigma32"    , "label": 1}),
            "sigma38": iPromoter2LDatabase(path.join(self.root_path    , "sigma38.fasta")    , {"tag": "sigma38"    , "label": 1}),
            "sigma54": iPromoter2LDatabase(path.join(self.root_path    , "sigma54.fasta")    , {"tag": "sigma54"    , "label": 1}),
            "sigma70": iPromoter2LDatabase(path.join(self.root_path    , "sigma70.fasta")    , {"tag": "sigma70"    , "label": 1}),
            "nonPromoter": iPromoter2LDatabase(path.join(self.root_path, "nonPromoter.fasta"), {"tag": "nonPromoter", "label": 0}),
        }
        self.chlg_model = chlg_model
        self.output_dir = output_dir
        os.makedirs( output_dir , exist_ok=True )
        self.output_file = "iPromoter2L-predictions.csv"
        self.pred_eval_fn = np.mean #max

    def run(self):
        # self.print_fn(self.databases)
        # self.generateMasterTable()
        # self.load_model(path.join("models", "{}.model".format(self.chlg_model)))
        # self.predict()
        self.evaluate()
        self.evaluate_per_sigma()
    
    def shuffle(self, df):
        return df.sample(frac=1).reset_index(drop=True)
        
    def generateMasterTable(self):
        tables = list()
        for i_d, key in enumerate(self.databases.keys()):
            database = self.databases[key]
            self.print_fn("BUILDING TABLE: {} AT {}".format(key, database.filepath))
            data = database.build()
            tables.append(data)
        tablesDF = pd.concat(tables)
        tableShuffled = self.shuffle(tablesDF)
        
        self.data = tableShuffled
        return tableShuffled

    def load_model(self, model_path):
        self.print_fn("\n\n LOADING ML MODEL {}".format(model_path)) 
        self.model = joblib.load(model_path)

    def window(self, input_seq, window_size=40):
        seqs = list()
        seq_len = len(input_seq)
        for i in range(seq_len - window_size + 1):
            seqs.append(input_seq[i:i+window_size].upper())
        return np.array(seqs)

    def predict(self):
        self.print_fn("\n\n CREATING PREDICTIONS USING MODEL {} ON DATA WITH SHAPE {}".format(
            self.chlg_model, self.data.shape))
        columns = list(self.data.columns) + ["predictions"]
        predictions = pd.DataFrame(columns=columns)
        for index, row in self.data.iterrows():
            self.print_fn("\n\n{}/{}. PREDICTING SEQ: {}".format(index+1, len(self.data), row["seq"])) 
            seqs = self.window(row["seq"])
            X = dataConverter(seqs, print_fn=self.silent_print_fn, data_type=self.chlg_model,
                              log_file=None, tokenizer_path="models/tokenizer.data")
            y_probs = self.model.predict_proba( X )
            y_pred  = y_probs[:, 1] if y_probs.shape[1] == 2 else y_probs[:, 0]
            dict_data = dict(row.copy())
            dict_data.update({"predictions": "-".join([str(x) for x in y_pred])})
            predictions = predictions.append(dict_data, ignore_index=True)
        predictions.to_csv(path.join(self.output_dir, self.output_file))
        self.predictions = predictions
        return predictions

    def evaluate(self):
        self.predictions = pd.read_csv(path.join(self.output_dir, self.output_file))
        y_pred = np.array([self.pred_eval_fn(np.array(pred.split("-")).astype(float)) for pred in self.predictions.predictions])
        y = self.predictions.label.values
        self.print_fn("Y_PRED ({}): {}".format(y_pred.shape, y_pred))
        self.print_fn("Y      ({}): {}".format(y.shape, y))
        ave_pre = average_precision_score(y, y_pred)
        roc_auc = roc_auc_score(y, y_pred)
        p, r, pr_t = precision_recall_curve(y, y_pred)
        fpr, tpr, roc_t = roc_curve(y, y_pred)

        self.print_fn("EVALUATION({}): \nAVE_PRE: {} ROC_AUC: {}".format("combined", ave_pre, roc_auc))
        data = {
            "p"   : list(p.astype(str)),
            "r"   : list(r.astype(str)),
            "pr_t": list(pr_t.astype(str)),
            "fpr" : list(fpr.astype(str)),
            "tpr" : list(tpr.astype(str)),
            "roc_t"  : list(roc_t.astype(str)),
            "ave_pre": ave_pre,
            "roc_auc": roc_auc
        }
        try:
            summary = pd.DataFrame({'p': p, 'r': r, 't': list(pr_t)+[1]})
            best_thres = summary[summary.p == summary.r].iloc[0].t
            mcc = matthews_corrcoef(y, y_pred >= best_thres)
            acc = accuracy_score(y, y_pred >= best_thres)
            f1  = f1_score(y, y_pred >= best_thres)

            data["best_threshold"] = {
                "t": summary[summary.p == summary.r].iloc[0].t,
                "p": summary[summary.p == summary.r].iloc[0].p,
                "r": summary[summary.p == summary.r].iloc[0].r,
                "acc": acc,
                "mcc": mcc,
                "f1" : f1
            }
            print("Best Threshold Value: ", data["best_threshold"])
        except Exception as e:
            print("ERROR: ", e)
        json.dump(data, open(path.join(self.output_dir, self.output_file.replace("csv", "json")), "w"), indent=2)

    def evaluate_per_sigma(self):
        self.predictions = pd.read_csv(path.join(self.output_dir, self.output_file))

        for i_s, sigma in enumerate(self.predictions.tag.unique()):
            notes = list()
            self.print_fn("{}/{}. EVALUATING SIGMA: {} VS NON-PROMOTERS".format(
                i_s+1, len(self.predictions.tag.unique()), sigma))
            data = self.predictions
            data = data[(data.tag == sigma) | (data.tag == self.NON_PROMOTERS_TAG)]
            self.print_fn("{}: P: {}/{} N: {}/{} TOTAL: {}/{}".format(
                sigma, len(data[data.label == 1]), len(data),
                len(data[data.label == 0]), len(data),
                len(data), len(self.predictions)))

            y_pred = np.array([self.pred_eval_fn(np.array(pred.split("-")).astype(float)) for pred in data.predictions])
            y = data.label.values
            self.print_fn("Y_PRED ({}): {}".format(y_pred.shape, y_pred))
            self.print_fn("Y      ({}): {}".format(y.shape, y))
            ave_pre = average_precision_score(y, y_pred)
            roc_auc = 0
            
            try:
                roc_auc = roc_auc_score(y, y_pred)
            except Exception as e: 
                notes.append(str(e))
            p, r, pr_t = precision_recall_curve(y, y_pred)
            fpr, tpr, roc_t = roc_curve(y, y_pred)

            print("EVALUATION({}): \nAVE_PRE: {} ROC_AUC: {}".format(sigma, ave_pre,  roc_auc))
            data = {
                "p"   : list(p.astype(str)),
                "r"   : list(r.astype(str)),
                "pr_t": list(pr_t.astype(str)),
                "fpr" : list(fpr.astype(str)),
                "tpr" : list(tpr.astype(str)),
                "roc_t"  : list(roc_t.astype(str)),
                "ave_pre": ave_pre,
                "roc_auc": roc_auc,
                "notes"  : str(notes)
            }
            try:
                summary = pd.DataFrame({'p': p, 'r': r, 't': list(pr_t) + [1]})
                best_thres = summary[summary.p == summary.r].iloc[0].t
                mcc = matthews_corrcoef(y, y_pred >= best_thres)
                acc = accuracy_score(y, y_pred >= best_thres)
                f1 = f1_score(y, y_pred >= best_thres)

                data["best_threshold"] = {
                    "t": summary[summary.p == summary.r].iloc[0].t,
                    "p": summary[summary.p == summary.r].iloc[0].p,
                    "r": summary[summary.p == summary.r].iloc[0].r,
                    "acc": acc,
                    "mcc": mcc,
                    "f1": f1
                }
                print("Best Threshold Value: ", data["best_threshold"])
            except Exception as e:
                print("ERROR: ", e)
            json.dump(data, open(path.join(self.output_dir, sigma+"_"+self.output_file.replace("csv", "json")), "w"), indent=2)

    def print_fn(self, msg, log_file=None, debug=True):
        if debug:
            print(msg)

    def silent_print_fn(self, msg, log_file=None):
        pass



