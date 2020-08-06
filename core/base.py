
import pandas as pd
import numpy  as np 
import os
import time
from pathlib import Path
import subprocess
import pickle
import sys
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['interactive'] == True
from pybedtools import BedTool
from Bio.Seq import Seq 
import progressbar
from utils import tetranucleotide_list_to_string_list, tetraToHotEncoding, promoterToTetraFreq
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from utils import h1

import tensorflow as tf
import joblib 


class MLModel():
  def __init__(self, model=None, cv=10, param_grid=[], title="Model", isFinal=False, existentModel=None ):
    self.title = title
    self.folds = cv
    if(existentModel != None):
      print("LOADING EXISTENT MODEL", existentModel)
    gridCVModel = existentModel if(existentModel != None) else  GridSearchCV(
        model, 
        param_grid, 
        cv      = cv, 
        iid     = False,
        scoring = ['average_precision','precision', 'recall'],
        refit   = 'average_precision',
        verbose = 2
    )
    self.model  = model if isFinal else gridCVModel
    print("CREATING ","FINAL" if isFinal else "NON-FINAL" , " ML-MODEL.NAME-", title ,"\n\n", self.model, "\n\n" )
    print("Model Parameters: ", param_grid)
  
  def train(self, X, y):
    print("TRAIN MODEL: ", self.title)
    self.history = self.model.fit(X,y)
    print("HISTORY\n{}\n\n".format(self.history))
  
  def drawCurves(self, X, y):
    print("DRAW CURVES: ", self.title)
    ave_precision = []
    ave_recall    = []
    ave_fpr       = []
    ave_tpr       = []
    ave_pre_score = []
    ave_roc_auc   = []

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    f.suptitle(self.title)
    man = plt.get_current_fig_manager()
    man.canvas.set_window_title(self.title)
    ax1.set_title('ROC PER FOLD')
    ax2.set_title('AUPRC PER FOLD')
    # Stratified KFold to distribute positives and negatives equally in each fold
    skf = StratifiedKFold(n_splits=self.folds )#n_split)
    skf.get_n_splits(X, y)
    # The algorithms sklearn roc and auprc curve algoritms output 
    # precision, recall, fpr, tpr with different shapes
    # To calculate the average of all we create a bi-dimensional matrix of each
    # With the size of the row of the minimum array output by the methods of sklearn
    min_pre_row_size = 999
    min_rec_row_size = 999
    min_fpr_row_size = 999
    min_tpr_row_size = 999
    fold_count = 0
    for train_index, test_index in skf.split(X, y):
      _, X_test = X[train_index], X[test_index]
      _, y_test = y[train_index], y[test_index]

      # Obtain Prediction & Performance Metrics
      predict_proba        = self.model.predict_proba(X_test, verbose=1)
      y_score              = predict_proba[:, 1] if predict_proba.shape[1] == 2 else predict_proba[:, 0]
      average_precision    = average_precision_score(y_test, y_score)
      precision, recall, _ = precision_recall_curve(y_test, y_score)
      fpr, tpr, _          = roc_curve(y_test,  y_score)
      roc_auc              = roc_auc_score(y_test,  y_score)
      # Save to memory
      ave_precision.append(precision)
      ave_recall.append(recall)
      ave_fpr.append(fpr)
      ave_tpr.append(tpr)
      ave_pre_score.append(average_precision)
      ave_roc_auc.append(roc_auc)
      # Get the minimum size of array per each metric
      if min_pre_row_size > len(precision):
        min_pre_row_size = len(precision)
      if min_rec_row_size > len(recall):
        min_rec_row_size = len(recall)
      if min_fpr_row_size > len(fpr):
        min_fpr_row_size = len(fpr)
      if min_tpr_row_size > len(tpr):
        min_tpr_row_size = len(tpr)
      fold_count = fold_count + 1 
    # Calculate AUPRC and AUROC
    np_ave_pre = np.array([]).reshape(0, min_pre_row_size)
    np_ave_rec = np.array([]).reshape(0, min_rec_row_size)
    np_ave_fpr = np.array([]).reshape(0, min_fpr_row_size)
    np_ave_tpr = np.array([]).reshape(0, min_tpr_row_size)
    for i in range(self.folds):
      ax1.plot(ave_fpr[i]   ,ave_tpr[i]       , alpha=0.4, label="F-{}-AUROC-{:.4f}".format(  i, ave_roc_auc[i]     ))
      ax2.plot(ave_recall[i], ave_precision[i], alpha=0.4, label="F-{}-AUPRC-{:.4f}".format(  i, ave_pre_score[i]   ))
      np_ave_pre = np.vstack( [np_ave_pre, ave_precision[i][:min_pre_row_size]] )
      np_ave_rec = np.vstack( [np_ave_rec, ave_recall[i][:min_rec_row_size]] )
      np_ave_fpr = np.vstack( [np_ave_fpr, ave_fpr[i][:min_fpr_row_size]] )
      np_ave_tpr = np.vstack( [np_ave_tpr, ave_tpr[i][:min_tpr_row_size]] )
    mean_tpr     = np_ave_tpr.mean(0)
    mean_fpr     = np_ave_fpr.mean(0) 
    ax1.plot(mean_fpr          , mean_tpr          , alpha=1.0 )
    ax2.plot(np_ave_rec.mean(0), np_ave_pre.mean(0), alpha=1.0 )
    ax1.plot([0, 1], [0, 1]    , linestyle='--', lw=2, color='r', label='Random', alpha=.8)
    ax2.plot([0, 1], [0.5, 0.5], linestyle='--', lw=2, color='r', label='Random', alpha=.8)
    leg1 = ax1.legend(loc='center right', bbox_to_anchor=(-0.2, 0.5))
    leg2 = ax2.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))
    f.suptitle("Params: "+ str(self.model.best_params_ ))
    ax1.set_title('AUROC: {:.4f}'.format(np.mean(ave_roc_auc))) 
    ax2.set_title('AUPRC: {:.4f}'.format(np.mean(ave_pre_score)))
    plt.draw()
    plt.pause(0.01)
    plt.pause(0.01)
    plt.savefig(self.title+'.png', bbox_extra_artists=[leg1, leg2] , bbox_inches='tight'  )
    return np.mean(ave_pre_score)
    
  def predict(self, X):
    predict_proba  = self.model.predict_proba(X)
    predict        = self.model.predict(X)
    return predict, predict_proba