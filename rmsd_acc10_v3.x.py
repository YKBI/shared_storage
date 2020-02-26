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
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.callbacks import EarlyStopping
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense, Activation

def help():
	print "print help usage\n"
	print "Usage: python %s [ peptide seq] [ pdb ]\n"%sys.argv[0]

if len(sys.argv) == 1:
	help()
else :
	GEAR = '/awork06-1/neoscan_gear'
	LIB = '/awork06-1/neoscan_lib'
	sam = sys.argv[1]
	sam1 = sys.argv[2]
	sam2 = sys.argv[3]

	try:
		if not os.path.exists('result') :
			os.mkdir('result',0777)
	except OSError:
		pass

	pdbs = []
	with open(sam1,'r') as f:
		lines = f.readlines()
		for line in lines:
			pdbs.append(line[:-1])

	model= tf.keras.models.load_model(LIB + '/rmsd_acc.h5')
	for pdb in pdbs:
		if os.path.exists(sam + '/' + pdb + '_energy_matrix/total_env_filter.txt') :
			os.remove(sam + '/' + pdb + '_energy_matrix/total_env_filter.txt')
		if os.path.exists(sam + '/' + pdb + '_energy_matrix/full_env.txt') :
			cmd = GEAR + '/ext_corr_full ' + LIB + '/total_phi_psi.txt ' + sam + '/' + pdb + '_energy_matrix/full_env.txt 5 > ' + sam + '/' + pdb + '_energy_matrix/total_env_filter.txt'
		else:
			cmd = GEAR + '/ext_corr_total ' + LIB + '/total_phi_psi.txt ' + sam + '/' + pdb + '_energy_matrix/total_env.txt 5 > ' + sam + '/' + pdb + '_energy_matrix/total_env_filter.txt'
		subprocess.call(cmd,shell=True)
		with open(sam + '/' + pdb + '_energy_matrix/total_env_filter.txt','r') as f:
			nf = f.readlines()
	
		if len(nf) > 1 :
			df = pd.read_csv(sam + '/' + pdb + '_energy_matrix/total_env_filter.txt',sep='\t')

			head = ['PDB','total_score','P1','P2','P3','P4','P5','P6','P7','P8','P9','PHI2','PHI3','PHI7','PSI2','PSI3','PSI7']
			feat = ['P1','P2','P3','P4','P5','P6','P7','P8','P9','PHI2','PHI3','PHI7','PSI2','PSI3','PSI7']
			head1 = ['PDB','total_score','P1','P2','P3','P4','P5','P6','P7','P8','P9','PHI2','PHI3','PHI7','PSI2','PSI3','PSI7','pred_rmsBB']
			df1 = df[head]
			X = df[feat]

	#		model= tf.keras.models.load_model(LIB + '/rmsd_acc.h5')
			predict_res = model.predict(X).flatten()
			df2 = pd.DataFrame()
			df2['pred_rmsBB']= predict_res
			df['pred_rmsBB'] = df2['pred_rmsBB']
			df[head1].to_csv('result/pred_' + sam + '_' + pdb + '_total_env.txt',sep='\t',index=False)
			cmd =  GEAR + '/ext_corr_dist ' + LIB + '/total_phi_psi.txt ' + 'result/pred_' + sam + '_' + pdb + '_total_env.txt ' + sam +  ' ' + sam2 + ' > result/pred_' + sam + '_' + pdb + '_total_env1.txt'
			subprocess.call(cmd,shell=True)
		else :
			print 'In %s_%s, No result in angle filrered result'%(sam,pdb)

#	plt.clf()
#	plt.figure()
#	X = df['rmsBB']
#	Y = df['pred_rmsBB']
#	X = sm.add_constant(X)
#	smodel = sm.OLS(Y, X).fit()
#	spredictions = smodel.predict(X)
#	text = 'test Fit,R-squared=%f'%(smodel.rsquared)
#	print text
#	plt.scatter(df['rmsBB'], df['pred_rmsBB'],label='original_data')
#	plt.plot(df['rmsBB'],spredictions,'r-',label=text)
#	plt.xlabel('True Values [rmsBB]')
#	plt.ylabel('Predictions [rmsBB]')
#	plt.legend()
#	plt.savefig(sam + '_test_fit.png')
