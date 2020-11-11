clc
%%%%%For the first layer, we selected the feature combination KNN(15)+BPB(130)+DNC(9)+MNC(1)+DAC(10)

disp("_______________________LOADING TRAINING FEATURES_______________________")
load('allfeatureL1.mat')
load('F-scorezhiL1.mat')

TrainFeatureVector = allfeatureL1;
TrainLabel         = yapp;

disp(["TRAIN SHAPE: (", size(TrainFeatureVector), ") \nSAMPLE: ", TrainFeatureVector(1, 1:5) ])
disp(["LABEL SHAPE: (", size(TrainLabel), ") \nSAMPLE: ", TrainLabel(1) ])

tic
disp("_______________________START TRAINING_______________________")
% model=svmtrain(TrainLabel,TrainFeatureVector,'-c 32 -g 0.01056 -w1 1 -w-1 1');
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 32 -g 0.01056 -w1 1 -w-1 1 -b 1');
timeElapsed = toc;
disp(["TIME ELAPSED (seconds): ", timeElapsed])

disp("_______________________SAVING MODEL TO MODEL.MAT_______________________")
save('model.mat', 'model')