import pandas as pd, numpy as np, time, progressbar
from pathlib import Path
import joblib, traceback, os, sys
from benchmark.ipromoter2l import iPromoter2L

def run_benchmark(model_tag, chlg_model, log_file=None, print_fn=print, output_dir="RF-HOT"):
    if log_file is None:
        print("NO LOG FILE SPECIFIED. REDIRECTING OUTPUT TO CONSOLE.")
        print_fn = print

    print_fn("RUNNING BENCHMARK ON MODEL: {}".format(model_tag))

    if model_tag == "iPromoter2L":
        model = iPromoter2L(chlg_model=chlg_model, output_dir=output_dir)
        model.run()

    
