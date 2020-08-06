import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from base_model import MLModel
# import keras
# from keras.models import Sequential
# from keras.layers import Dense, Embedding, LSTM, GRU
# from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.utils import class_weight
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Embedding, LSTM, GRU, Activation, Dropout
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from tensorflow.keras.callbacks import EarlyStopping,  ModelCheckpoint, CSVLogger
from tensorflow.keras.metrics import TruePositives, FalsePositives, TrueNegatives, FalseNegatives, BinaryAccuracy, Precision, Recall, AUC
import json
import copy
from utils import testTry
import joblib

class RandomForest(MLModel):
  def __init__(self, class_weights, param_grid):
    self.param_grid = param_grid
    print("CLASS WEIGHTS: ", dict(enumerate(class_weights)))
    MLModel.__init__(
        self, 
        model=RandomForestClassifier(
          min_samples_leaf=5, 
          # n_estimators=500
          verbose=2,
          class_weight=dict(enumerate(class_weights))
        ), 
        param_grid=self.param_grid, 
        title="Random Forest")
  def train(self, X, y):
    print("TRAIN RF MODEL: ", self.title)
    class_weights = class_weight.compute_class_weight('balanced', np.unique(y.flatten()), y.flatten())
    print("""
    Class weights: \n{}\n{}\n
    for classes: \n{}\n
    # Promoters: {}
    # of Non-Promoters: {}
    % Promoters: {}%
    % of Non-Promoters: {}%
    """.format( 
        class_weights, 
        dict(enumerate(class_weights)),
        np.unique(y.flatten()), 
        np.count_nonzero(y == 1), 
        np.count_nonzero(y == 0),
        np.count_nonzero(y == 1)*100/len(y), 
        np.count_nonzero(y == 0)*100/len(y),
    ))
    print("""INPUTS
    X:      {}
    y:      {}
    """.format(X.shape, y.shape ))
    self.history = self.model.fit(
      X,y, 
      sample_weight   = class_weights,
    )
    print("\n\n{} {} {}\n\n".format( 10*"_ " , "RF TRAINING RESULTS" , 10*"_ "))
    print("""
      BEST ESTIMATOR:          {} 
      BEST SCORE:              {}
      BEST PARAMS:             {}
      BEST INDEX IN CV SEARCH: {}
      SCORER FUNCTIONS:        {}
      \n
      HISTORY OBJ:             {}        
    \n\n""".format( 
      self.history.best_estimator_,
      self.history.best_score_ , 
      self.history.best_params_ ,
      self.history.best_index_,
      self.history.scorer_,
      self.history
    ))
    print("cv_results_dict: ")
    print(pd.DataFrame( self.history.cv_results_ ))
    print("\n\nSAVING HISTORY: {}.history\n\n".format(self.title))
    try:
      historyDict = copy.copy(self.history.__dict__)
      print(historyDict.keys())
      testTry(lambda : historyDict.pop("best_estimator_") , "" )
      testTry(lambda : historyDict.pop("multimetric_") , "" )
      testTry(lambda : historyDict.pop("scorer_") , "" )
      testTry(lambda : historyDict.pop("estimator")    , "" )
      testTry(lambda : historyDict["best_params_"].pop('keras_eval_metric') , "best_params_" )
      testTry(lambda : historyDict["cv_results_"].pop('param_dropout'), "cv_results_" )
      testTry(lambda : historyDict["cv_results_"].pop('param_dropout_rate'), "cv_results_" )
      testTry(lambda : historyDict["cv_results_"].pop('param_keras_eval_metric'), "cv_results_" )
      testTry(lambda : historyDict["param_grid"].pop('keras_eval_metric'), "param_grid" )
      joblib.dump(  historyDict, "{}.history".format(self.title) )
      print("\n\nHISTORY: {}.history SAVED SUCCESSFULLY.\n\n".format(self.title))
    except Exception as e:
      print("\n\nCANNOT SAVE {} BECAUSE: {}\n\n".format(self.title , str(e)))
    return self.history

def create_model_RNN_GRU( 
  optimizer='adam', units=200, activation='sigmoid', 
  EMBEDDING_DIM=25, max_length=37,  vocab_size=350, hidden_dims=1, hidden_dim_units=25 ):
  keras_eval_metric = [[ 
    TruePositives(name='tp'),
    FalsePositives(name='fp'),
    TrueNegatives(name='tn'),
    FalseNegatives(name='fn'), 
    BinaryAccuracy(name='accuracy'),
    Precision(name='precision'),
    Recall(name='recall'),
    AUC(name='auc'),
  ]]
  model = Sequential()
  print("Parameters", "units", units, "activation", activation, "EMBEDDING_DIM", EMBEDDING_DIM, "max_length", max_length,  "vocab_size", vocab_size, " hidden_dims  ", hidden_dims )
  model.add(Embedding(vocab_size, EMBEDDING_DIM,  input_length=max_length))
  model.add(GRU(units=hidden_dim_units, dropout=0.2, recurrent_dropout=0.2, activation=activation ))
  for i in range(hidden_dims):
    model.add(Dense(hidden_dim_units))
    model.add(Activation(activation))
    model.add(Dropout(0.2))
  model.add(Dense(1))
  model.add(Activation(activation))
  model.compile(loss='binary_crossentropy', optimizer=optimizer , metrics=keras_eval_metric  )
  print( "\n\n", model.summary(), "\n\n", model.get_config(), "\n\n")
  return model


def create_model_RNN_LSTM( 
    optimizer='adam', units=200, activation='sigmoid', 
    EMBEDDING_DIM=25, max_length=37,  vocab_size=350, hidden_dims=1, hidden_dim_units=25 ):
    print("Parameters", "units", units, "activation", activation, "EMBEDDING_DIM", EMBEDDING_DIM, "max_length", max_length,  "vocab_size", vocab_size, " hidden_dims  ", hidden_dims )
    keras_eval_metric = [[ 
      TruePositives(name='tp'),
      FalsePositives(name='fp'),
      TrueNegatives(name='tn'),
      FalseNegatives(name='fn'), 
      BinaryAccuracy(name='accuracy'),
      Precision(name='precision'),
      Recall(name='recall'),
      AUC(name='auc'),
    ]]
    model = Sequential()
    model.add(Embedding(vocab_size, EMBEDDING_DIM,  input_length=max_length))
    model.add(LSTM(units=hidden_dim_units, dropout=0.2, recurrent_dropout=0.2, activation=activation ))
    for i in range(hidden_dims):
      model.add(Dense(hidden_dim_units))
      model.add(Activation(activation))
      model.add(Dropout(0.2))
    model.add(Dense(1))
    model.add(Activation(activation))
    model.compile(loss='binary_crossentropy', optimizer=optimizer , metrics=keras_eval_metric  )
    print( "\n\n", model.summary(), "\n\n", model.get_config(), "\n\n")

    return model

def getRNNGrid():
  seed = 7
  np.random.seed(seed)
  hidden_dim_units=[100]
  EMBEDDING_DIM=[50]
  param_grid = dict(hidden_dim_units=hidden_dim_units, EMBEDDING_DIM=EMBEDDING_DIM)  
  print("\nCREATING PARAM GRID: \n\n{}\n\n".format(param_grid))
  #Other params
  #optimizer    = ['SGD', 'RMSprop', 'Adagrad', 'Adadelta', 'Adam', 'Adamax', 'Nadam']
  #hidden_layer = [128, 256, 512, 1024, 2048]
  #activation   = ["softmax", 'sigmoid']
  #param_grid   = dict(optimizer=optimizer, hidden_layer=hidden_layer, activation=activation)  
  return param_grid

class RNN(MLModel):
  def __init__(self, vocab_size, create_model_callback, title, layers, isFinal, existentModel ):
    print("""\n\nRNN PARAMETERS
    _________________________________
    vocab_size:    {}
    title:         {}
    layers:        {}
    isFinal:       {} 
    existentModel: {} 
    \n\n""".format(
      vocab_size, 
      title, 
      layers, 
      isFinal, 
      existentModel
    ))
    metric_monitor = "val_recall" #"val_loss"
    self.callbacks   =  [ 
      EarlyStopping(monitor=metric_monitor, mode='min', verbose=1) , 
      ModelCheckpoint( "{}_checkpoint_model.h5".format(title),monitor=metric_monitor , mode='max', save_best_only=True, verbose=1),
      CSVLogger('{}_train_callback_log.txt'.format(title))
    ]
    MLModel.__init__(
        self, 
        model= KerasClassifier(
          build_fn=create_model_callback, epochs=50, batch_size=10, 
          verbose=2, vocab_size=vocab_size, hidden_dims=layers,
        ), 
        param_grid=getRNNGrid(), 
        title=title, 
        isFinal=isFinal,
        existentModel=existentModel)

  def train(self, X, y, X_test, y_test):
    print("TRAIN RNN MODEL: ", self.title)
    class_weights = class_weight.compute_class_weight('balanced', np.unique(y.flatten()), y.flatten())
    print("""
    Class weights: \n{}\n{}\n
    for classes: \n{}\n
    # Promoters: {}
    # of Non-Promoters: {}
    % Promoters: {}%
    % of Non-Promoters: {}%
    """.format( 
        class_weights, 
        dict(enumerate(class_weights)),
        np.unique(y.flatten()), 
        np.count_nonzero(y == 1), 
        np.count_nonzero(y == 0),
        np.count_nonzero(y == 1)*100/len(y), 
        np.count_nonzero(y == 0)*100/len(y),
    ))
    print("""INPUTS
    X:      {}
    y:      {}
    X_test: {}
    y_test: {}
    """.format(X.shape, y.shape, X_test.shape, y_test.shape ))
    self.history = self.model.fit(
      X,y, 
      callbacks= self.callbacks,
      validation_data= (X_test, y_test),
      class_weight   = class_weights,
    )
    print("\n\n{} {} {}\n\n".format( 10*"_ " , "RNN TRAINING RESULTS" , 10*"_ "))
    print("""
      BEST ESTIMATOR:          {} 
      BEST SCORE:              {}
      BEST PARAMS:             {}
      BEST INDEX IN CV SEARCH: {}
      SCORER FUNCTIONS:        {}
      \n
      HISTORY OBJ:             {}        
    \n\n""".format( 
      self.history.best_estimator_,
      self.history.best_score_ , 
      self.history.best_params_ ,
      self.history.best_index_,
      self.history.scorer_,
      self.history
    ))

    try:
      print("SAVING RNN MODEL BEST ESTIMATOR : ", self.title)
      self.model.best_estimator_.model.save( "{}_best_estimator.h5".format(self.title) )
    except Exception as e:
      print("Couldn't save best estimator: ", str(e))

    print("cv_results_dict: ")
    print(pd.DataFrame( self.history.cv_results_ ))
    print("\n\nSAVING HISTORY: {}.history\n\n".format(self.title))
    try:
      historyDict = copy.copy(self.history.__dict__)
      print(historyDict.keys())
      testTry(lambda : historyDict.pop("best_estimator_") , "" )
      testTry(lambda : historyDict.pop("multimetric_") , "" )
      testTry(lambda : historyDict.pop("scorer_") , "" )
      testTry(lambda : historyDict.pop("estimator")    , "" )
      testTry(lambda : historyDict["best_params_"].pop('keras_eval_metric') , "best_params_" )
      testTry(lambda : historyDict["cv_results_"].pop('param_dropout'), "cv_results_" )
      testTry(lambda : historyDict["cv_results_"].pop('param_dropout_rate'), "cv_results_" )
      testTry(lambda : historyDict["cv_results_"].pop('param_keras_eval_metric'), "cv_results_" )
      testTry(lambda : historyDict["param_grid"].pop('keras_eval_metric'), "param_grid" )
      joblib.dump(  historyDict, "{}.history".format(self.title) )
      print("\n\nHISTORY: {}.history SAVED SUCCESSFULLY.\n\n".format(self.title))
    except Exception as e:
      print("\n\nCANNOT SAVE {} BECAUSE: {}\n\n".format(self.title , str(e)))
    return self.history