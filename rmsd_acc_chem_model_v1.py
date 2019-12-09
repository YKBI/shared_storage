#!/usr/bin/env python

import os
import sys
import subprocess
import shutil
import string
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
import statsmodels.api as sm
import sklearn
from numpy.random import seed
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from tensorflow.python.keras.callbacks import TensorBoard
from keras.callbacks import EarlyStopping
from keras.models import Sequential
from keras.models import load_model
from keras import models, layers
from keras.layers import Dense, Dropout, Activation
from time import time

#seed(10)

def fit_figure(sam,real,pred):
	plt.clf()
	plt.figure()
	X = real
	Y = pred
#	X = sm.add_constant(X.to_numpy())
	X = sm.add_constant(X)
	smodel = sm.OLS(Y, X).fit()
	spredictions = smodel.predict(X)
	res = Y - spredictions
	ss_tot = np.sum( (Y - Y.mean())**2 )
	ss_res = np.sum( (Y - spredictions)**2 )
	man_rsq = 1 - ss_res/ss_tot
	auto_rsq = sklearn.metrics.r2_score(Y, spredictions)
	lig = sam.split('_')[0]
	text = '%s Fit,R-squared=%f'%(lig,smodel.rsquared)
	text1 = '%s Fit,manually calc.R-squared=%f'%(lig,man_rsq)
	text2 = '%s Fit,sklearn R-squared=%f'%(lig,auto_rsq)
	with open('total_rsq.txt','a') as rf:
		rf.write('%s\n'%(text))
		rf.write('%s\n'%(text1))
		rf.write('%s\n'%(text2))
	print text
	print text1
	print text2
	plt.scatter(real, pred,label='original_data')
	plt.plot(real,spredictions,'r-',label=text)
	plt.xlabel('True Values [rmsBB]')
	plt.ylabel('Predictions [rmsBB]')
	plt.legend()
	plt.savefig(sam + '_test_fit.png')

def plot_history(history):
	hist = pd.DataFrame(history.history)
	hist['epoch'] = history.epoch

	plt.figure(figsize=(8,12))
	
	plt.subplot(2,1,1)
	plt.xlabel('Epoch')
	plt.ylabel('Mean Abs Error')
	plt.plot(hist['epoch'], hist['mean_absolute_error'],label='Train Error')
	plt.plot(hist['epoch'], hist['val_mean_absolute_error'],label='Val Error')
	plt.legend()

	plt.subplot(2,1,2)
	plt.xlabel('Epoch')
	plt.ylabel('Mean Square Error [rmsBB^2]')
	plt.plot(hist['epoch'], hist['mean_squared_error'],label='Train Error')
	plt.plot(hist['epoch'], hist['val_mean_squared_error'],label='Val Error')
	plt.legend()

	plt.savefig(sam + '_qual.png')

sam = sys.argv[1]
sam1 = sys.argv[2]

df = pd.read_csv(sam + '_full_rec_env_dist_1.txt',sep='\t')

if len(sam.split('_')) > 2 :
	lig = sam.split('_')[4]
	rec = sam.split('_')[0]
else:
	lig = sam.split('_')[1]
	rec = sam.split('_')[0]

feats = []
feats.append('ligand_rms_no_super_X')
with open(sam + '_feat','r') as fx:
	lines = fx.readlines()
	for line in lines:
		feats.append(line[:-1])
#feats = df.columns.tolist()

#feats.remove('PDB')
#feats.remove('rmsd')
df1 = df[feats]
#df1.to_csv('test.txt',sep='\t',index=False)
train_df = df1.sample(frac=0.8,random_state=0)
test_df = df1.drop(train_df.index)


train_labels =  train_df.pop('ligand_rms_no_super_X')
test_labels = test_df.pop('ligand_rms_no_super_X')

#lr = LinearRegression().fit(train_df,train_labels)
#ridge = Ridge(alpha=10).fit(train_df,train_labels)

#print 'Training score :%f'%(ridge.score(train_df,train_labels))
#print 'Test score: %f'%(ridge.score(test_df,test_labels))

def build_model():
	model = Sequential()
	model.add(layers.Dense(64,input_dim=len(feats)-1,activation='relu'))
	for i in range(int(sam1)):
		model.add(layers.Dense(64,activation='relu'))
	model.add(layers.Dense(1,activation='relu'))
	model.compile(loss='mean_squared_error',optimizer='adam',metrics=['mean_absolute_error', 'mean_squared_error'])
	return model

model = build_model()

model.summary()

early_stopping = EarlyStopping(patience=20)
history = model.fit(train_df, train_labels, epochs=1000, batch_size=128,validation_split = 0.2, callbacks=[early_stopping])

plot_history(history)

loss_and_metrics = model.evaluate(test_df,test_labels)

print loss_and_metrics

model.save(rec + '_' + lig + '_model.h5')

pred_res = model.predict(test_df).flatten()

fit_figure(sam,test_labels,pred_res)
