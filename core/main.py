from data import generateData, loadData, saveData
from models import RandomForest, RNN, create_model_RNN_GRU, create_model_RNN_LSTM
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import argparse
import numpy as np
import os
import time
from keras.wrappers.scikit_learn import KerasClassifier
from base_model import MLModel
# from sklearn.externals import joblib
 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
tf.compat.v1.logging.set_verbosity( tf.compat.v1.logging.ERROR)
import joblib
from sklearn.utils import class_weight

train_unbalanced = True

def trainWithHotEncoding(
  hot_encoded_train_features, hot_encoded_train_labels, 
  hot_encoded_test_features, hot_encoded_test_labels,
  results, algorithms, isTesting
):
  # TRAIIN RANDOM FOREST
  class_weights = class_weight.compute_class_weight(
    'balanced', 
    np.unique(hot_encoded_train_labels.flatten()), 
    hot_encoded_train_labels.flatten()
  )
  param_grid = {  
    'max_features': [
      # None, 
      # "sqrt", 
      "log2" 
    ], 
    'n_estimators' : [
      # 1000, 
      2000, 
      # 3000
    ] 
  }
  if(isTesting):
    print("TESTING MODE RF: ONLY TRAINING 1 MODEL")
    param_grid={}

  rf_hot_encoding = RandomForest(class_weights=class_weights, param_grid=param_grid)
  rf_hot_encoding.title = "RANDOM FOREST HOT ENCODING TRAIN"
  rf_hot_encoding.train( hot_encoded_train_features, hot_encoded_train_labels)

  model_name = 'RF-HOT.model'
  print("SAVING MODEL: ", model_name )
  try:
    joblib.dump(rf_hot_encoding.model, model_name )
  except Exception as e:
    print("Cannot save {} because: \n\n".format(model_name) , str(e))

  algorithms["RANDOM FOREST HOT ENCODING"] = rf_hot_encoding
  rf_hot_encoding.drawCurves( X=hot_encoded_train_features,  y=hot_encoded_train_labels )
  rf_hot_encoding.title = "RANDOM FOREST HOT ENCODING TEST"
  results["RANDOM FOREST HOT ENCODING"] = rf_hot_encoding.drawCurves(
      X=hot_encoded_test_features, 
      y=hot_encoded_test_labels
  )
  return rf_hot_encoding, results, algorithms

def trainWithFrecuencies(
  tetra_freq_train_features, tetra_freq_train_labels, 
  tetra_freq_test_features, tetra_freq_test_labels, 
  results, algorithms, isTesting
):
  class_weights = class_weight.compute_class_weight(
    'balanced', 
    np.unique(tetra_freq_train_labels.flatten()), 
    tetra_freq_train_labels.flatten()
  )
  param_grid = {  
    'max_features': [
      # None, 
      # "sqrt", 
      "log2" 
    ], 
    'n_estimators' : [
      # 1000, 
      2000, 
      # 3000
    ] 
  }
  if(isTesting):
    print("TESTING MODE RF: ONLY TRAINING 1 MODEL")
    param_grid={}

  rf = RandomForest(class_weights=class_weights, param_grid=param_grid)
  rf.title = "RANDOM FOREST TETRA NUCLEOTIDE FREQUENCY TRAIN"
  rf.train(tetra_freq_train_features, tetra_freq_train_labels)
  
  model_name = 'RF-TETRA.model'
  print("SAVING MODEL USING JOBLIB: ", model_name )
  try:
    joblib.dump(rf.model, model_name )
  except Exception as e:
    print("Cannot save {} because: \n\n".format(model_name) , str(e))

  algorithms["RANDOM FOREST"] = rf
  rf.drawCurves( X=tetra_freq_train_features, y=tetra_freq_train_labels)
  rf.title = "RANDOM FOREST TETRA NUCLEOTIDE FREQUENCY TEST"
  results["RANDOM FOREST"] = rf.drawCurves( X=tetra_freq_test_features, y=tetra_freq_test_labels)
  return rf, results, algorithms

def trainWithRNNGRU(
  rnn_token_train_features, rnn_token_train_labels,  
  rnn_token_test_features, rnn_token_test_labels,
  results, algorithms, vocab_size, layers, isFinal, doTraining, existentModel
):
  print("T W GRU vocab", vocab_size)
  rnn      = RNN(vocab_size=vocab_size, create_model_callback=create_model_RNN_GRU, title="GRU-TRAIN-D-{}".format(layers) , layers=layers, isFinal=isFinal, existentModel=existentModel)
  print(
      np.array(rnn_token_train_features).shape, 
      np.array(rnn_token_train_labels).shape )

  if(doTraining):
    rnn.train(
        X=np.array(rnn_token_train_features, dtype=float), 
        y=np.array(rnn_token_train_labels  , dtype=float),
        X_test=np.array(rnn_token_test_features, dtype=float), 
        y_test=np.array(rnn_token_test_labels  , dtype=float)
    )
  else:
    print(rnn.title, "-NOT TRAINING")
    
  model_name = 'GRU-D-{}.model'.format(layers)
  print("SAVING MODEL & WEIGHTS USING H5: ", model_name )
  try:
    print(rnn.model.best_estimator_.model)
    rnn.model.best_estimator_.model.save_weights(model_name+"_weights.h5")
    rnn.model.best_estimator_.model.save(model_name+".h5")
  except Exception as e:
    print("Cannot save {} because: {}\n\n".format(model_name, str(e)))

  algorithms["RNN_GRU"] = rnn
  rnn.drawCurves(
    X=np.array(rnn_token_train_features, dtype=float), 
    y=np.array(rnn_token_train_labels  , dtype=float)
  )
  rnn.title = "GRU-TEST-D-{}".format(layers)
  results["RNN_GRU"] = rnn.drawCurves(
    X=np.array(rnn_token_test_features, dtype=float), 
    y=np.array(rnn_token_test_labels  , dtype=float)
  )

  return rnn , results, algorithms

def trainWithRNNLSTM(
  rnn_token_train_features, rnn_token_train_labels,  
  rnn_token_test_features, rnn_token_test_labels,
  results, algorithms, vocab_size, layers, isFinal, doTraining, existentModel
):
  rnn      = RNN(vocab_size=vocab_size, create_model_callback=create_model_RNN_LSTM, title="LSTM-TRAIN-D-{}".format(layers), layers=layers, isFinal=isFinal, existentModel=existentModel )
  print(
      np.array(rnn_token_train_features).shape, 
      np.array(rnn_token_train_labels).shape,
      "LSTM VOCAB", vocab_size )

  if(doTraining):
    rnn.train(
        X=np.array(rnn_token_train_features, dtype=float), 
        y=np.array(rnn_token_train_labels  , dtype=float),
        X_test=np.array(rnn_token_test_features, dtype=float), 
        y_test=np.array(rnn_token_test_labels  , dtype=float)
    )
  else:
    print(rnn.title, "-NOT TRAINING")

  model_name = 'LSTM-D-{}.model'.format(layers)
  print("SAVING MODEL & WEIGHTS USING H5: ", model_name )
  try:
    print(rnn.model.best_estimator_.model)
    rnn.model.best_estimator_.model.save_weights(model_name+".h5")
    rnn.model.best_estimator_.model.save(model_name+".h5")
  except Exception as e:
    print("Cannot save {} because: {}\n\n".format(model_name, str(e)))

  algorithms["RNN_LSTM"] = rnn
  rnn.drawCurves(
    X=np.array(rnn_token_train_features, dtype=float), 
    y=np.array(rnn_token_train_labels  , dtype=float)
  )
  rnn.title = "LSTM-TEST-D-{}".format(layers)
  results["RNN_LSTM"] = rnn.drawCurves(
    X=np.array(rnn_token_test_features, dtype=float), 
    y=np.array(rnn_token_test_labels  , dtype=float)
  )

  return rnn , results, algorithms

def main():
  results       = {}
  algorithms    = {}

  switcher = {
    "1": "Training Using Random Forest with Hot Encoded Features",
    "2": "Training Using Random Forest with Tetranucleotide Features",
    "3": "Training Using RNN GRU ",
    "4": "Training Using RNN LSTM",
  }
  parser = argparse.ArgumentParser(description='Train to Find Promoters')
  parser.add_argument('--algorithm')
  parser.add_argument('--generate')
  parser.add_argument('--layers')
  parser.add_argument("-t", '--testing', help="Reduce Data For Testing?", action="store_true")
  args = parser.parse_args()
  print(
    """\n\nTERMINAL PARAMETERS
    _________________________________
    ALGORITHM:      {} {}
    GENERATE DATA:  {}
    DENSE-LAYERS:   {}
    REDUCE DATA 
    FOR TESTING:    {}
    \n\n
    """.format(
    args.algorithm, switcher.get(args.algorithm, "Invalid Algorithm"),
    args.generate, 
    args.layers,
    args.testing
  ))
  
  if((args.algorithm=="3" or args.algorithm=="4") and (args.layers==None) ):
    print("No title or layer counts for RNN")
    return

  data_path= "./data/UNBALANCED/" if train_unbalanced else "./data/BALANCED/"
  print("TRAINING UNBALANCED? {} DATA PATH: {}".format(train_unbalanced, data_path))
  hot_encoded_train_features, \
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
  vocab_size, tokenizer = generateData(
    None
    # pd.DataFrame([
    #     ["ECOLI_2"],
    # ], columns=["BACTERIA"])
  ) if args.generate else loadData(data_path=data_path)

  if(args.testing):
    print(
    """\n\nORIGINAL DATA
    ___________________________
    hot_encoded_train_features   {}
    hot_encoded_train_labels     {}
    hot_encoded_test_features    {}
    hot_encoded_test_labels      {}
    tetra_freq_train_features       {}
    tetra_freq_train_labels         {}
    tetra_freq_test_features        {}
    tetra_freq_test_labels          {}
    rnn_token_train_features   {}
    rnn_token_train_labels     {}
    rnn_token_test_features    {}
    rnn_token_test_labels      {}
    vocab_size                  {}
    tokenizer                   {}
    Results                     {}
    Algorithms                  {}
    # Promoters                 {}
    # Non-promoters             {}           
    % Promoters                 {}%
    % Non-promoters             {}%
    \n\n""".format(
      hot_encoded_train_features.shape,
      hot_encoded_train_labels.shape,
      hot_encoded_test_features.shape,
      hot_encoded_test_labels.shape,
      tetra_freq_train_features.shape,
      tetra_freq_train_labels.shape,
      tetra_freq_test_features.shape,
      tetra_freq_test_labels.shape,
      rnn_token_train_features.shape,
      rnn_token_train_labels.shape,
      rnn_token_test_features.shape,
      rnn_token_test_labels.shape,
      vocab_size, 
      tokenizer,
      results, 
      algorithms ,
      len(hot_encoded_train_labels[hot_encoded_train_labels==1]),
      len(hot_encoded_train_labels[hot_encoded_train_labels==0]),
      len(hot_encoded_train_labels[hot_encoded_train_labels==1])*100/len(hot_encoded_train_labels),
      len(hot_encoded_train_labels[hot_encoded_train_labels==0])*100/len(hot_encoded_train_labels),
    ))
    reduced_data_size = 1000
    print("REDUCING DATA FOR TESTING TO SIZE OF: ", reduced_data_size)
    hot_encoded_train_features = hot_encoded_train_features[:reduced_data_size]
    hot_encoded_train_labels   = hot_encoded_train_labels[:reduced_data_size]
    hot_encoded_test_features  = hot_encoded_test_features[:reduced_data_size]
    hot_encoded_test_labels    = hot_encoded_test_labels[:reduced_data_size]
    tetra_freq_train_features     = tetra_freq_train_features[:reduced_data_size]
    tetra_freq_train_labels       = tetra_freq_train_labels[:reduced_data_size]
    tetra_freq_test_features      = tetra_freq_test_features[:reduced_data_size]
    tetra_freq_test_labels        = tetra_freq_test_labels[:reduced_data_size]
    rnn_token_train_features = rnn_token_train_features[:reduced_data_size]
    rnn_token_train_labels   = rnn_token_train_labels[:reduced_data_size]
    rnn_token_test_features  = rnn_token_test_features[:reduced_data_size]
    rnn_token_test_labels    = rnn_token_test_labels[:reduced_data_size]
  
  # vocab_size = 550
  start_time = time.time()
  print(
    """\n\nDATA & SHAPES
    ___________________________
    hot_encoded_train_features   {}
    hot_encoded_train_labels     {}
    hot_encoded_test_features    {}
    hot_encoded_test_labels      {}
    tetra_freq_train_features       {}
    tetra_freq_train_labels         {}
    tetra_freq_test_features        {}
    tetra_freq_test_labels          {}
    rnn_token_train_features   {}
    rnn_token_train_labels     {}
    rnn_token_test_features    {}
    rnn_token_test_labels      {}
    vocab_size                  {}
    tokenizer                   {}
    start_time                  {}
    Results                     {}
    Algorithms                  {}
    # Promoters                 {}
    # Non-promoters             {}           
    % Promoters                 {}%
    % Non-promoters             {}%
    \n\n""".format(
      hot_encoded_train_features.shape,
      hot_encoded_train_labels.shape,
      hot_encoded_test_features.shape,
      hot_encoded_test_labels.shape,
      tetra_freq_train_features.shape,
      tetra_freq_train_labels.shape,
      tetra_freq_test_features.shape,
      tetra_freq_test_labels.shape,
      rnn_token_train_features.shape,
      rnn_token_train_labels.shape,
      rnn_token_test_features.shape,
      rnn_token_test_labels.shape,
      vocab_size, 
      tokenizer,
      start_time ,
      results, 
      algorithms,
      len(hot_encoded_train_labels[hot_encoded_train_labels==1]),
      len(hot_encoded_train_labels[hot_encoded_train_labels==0]),
      len(hot_encoded_train_labels[hot_encoded_train_labels==1])*100/len(hot_encoded_train_labels),
      len(hot_encoded_train_labels[hot_encoded_train_labels==0])*100/len(hot_encoded_train_labels),
  ))


  if(args.algorithm and args.algorithm == "1"):
    rf_hot_encoding, results, algorithms = trainWithHotEncoding(
      hot_encoded_train_features, hot_encoded_train_labels, 
      hot_encoded_test_features, hot_encoded_test_labels,
      results, algorithms, args.testing
    )

  elif(args.algorithm and args.algorithm == "2"):
    rf, results, algorithms = trainWithFrecuencies(
      tetra_freq_train_features, tetra_freq_train_labels, 
      tetra_freq_test_features, tetra_freq_test_labels, 
      results, algorithms, args.testing
    )

  elif(args.algorithm and args.algorithm == "3"):
    rnn , results, algorithms = trainWithRNNGRU(
      rnn_token_train_features, rnn_token_train_labels,  
      rnn_token_test_features, rnn_token_test_labels,
      results, algorithms, vocab_size, int(args.layers), False, True, None
    )
  
  elif(args.algorithm and args.algorithm == "4"):
    rnnLSTM , resultsLSTM, algorithms = trainWithRNNLSTM(
      rnn_token_train_features, rnn_token_train_labels,  
      rnn_token_test_features, rnn_token_test_labels,
      results, algorithms, vocab_size, int(args.layers), False, True, None
    )

  elapsed_time = time.time() - start_time
  print("\n\nELAPSED TIME: {} sec. or {} min.".format(elapsed_time, elapsed_time/60) )

  plt.show()
  
  
if __name__ == "__main__":
  main()