#!/usr/bin/env python

import os
import sys
import subprocess
import shutil
import string
import glob
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
	print "Usage: python %s [ peptide list] [ pdb ] [ HLA-type] [ biomarker / Target_gene / TP53 or KRAS]\n"%sys.argv[0]

if len(sys.argv) == 1:
	help()
else :
	GEAR = '/awork06-1/neoscan_gear'
	LIB = '/awork06-1/neoscan_lib'
	sam = sys.argv[1]
	sam1 = sys.argv[2]
	sam2 = sys.argv[3]
	sam3 = sys.argv[4]

	try:
		if not os.path.exists('result') :
			os.mkdir('result',0777)
	except OSError:
		pass

	peps = []
	with open(sam,'r') as f1:
		lines = f1.readlines()
		for line in lines:
			peps.append(line[:-1])

	pdbs = []
	with open(sam1,'r') as f:
		lines = f.readlines()
		for line in lines:
			pdbs.append(line[:-1])

	head = ['PDB','total_score','P1','P2','P3','P4','P5','P6','P7','P8','P9','PHI2','PHI3','PHI7','PSI2','PSI3','PSI7']
	feat = ['P1','P2','P3','P4','P5','P6','P7','P8','P9','PHI2','PHI3','PHI7','PSI2','PSI3','PSI7']
	head1 = ['PDB','total_score','P1','P2','P3','P4','P5','P6','P7','P8','P9','PHI2','PHI3','PHI7','PSI2','PSI3','PSI7','pred_rmsBB']

	model= tf.keras.models.load_model(LIB + '/rmsd_acc.h5')

	wdir = os.getcwd()

	with open('result/pred_' + sam + '_' + sam3 + '_pre.txt','a') as rf :
		rf.write('Sequence\tHLA-type\tPDB\ttotal_score\tP1\tP2\tP3\tP4\tP5\tP6\tP7\tP8\tP9\tPHI2\tPHI3\tPHI7\tPSI2\tPSI3\tPSI7\tpred_rmsBB\tDist\tClass\n')

	print wdir
	# Calculation of pred_rmsd
	for pep in peps :

		try:
			if not os.path.exists(pep + '_res'):
				os.mkdir(pep + '_res',0777)
		except OSError:
			pass

		for pdb in pdbs :
			df = pd.read_csv(pep + '/' + pdb + '_energy_matrix/full_env.txt',sep='\t')

			df1 = df[head]
			X = df[feat]

			predict_res = model.predict(X).flatten()
			df2 = pd.DataFrame()
			df2['pred_rmsBB']= predict_res
			df['pred_rmsBB'] = df2['pred_rmsBB']
			df[head1].to_csv(pep + '_res/pred_' + pep + '_' + pdb + '_total_env.txt',sep='\t',index=False)
			# Calculation of angle rms difference
			cmd = GEAR + '/ext_corr_dist ' + LIB + '/total_phi_psi.txt ' + pep + '_res/pred_' + pep + '_' + pdb + '_total_env.txt ' + pep +  ' ' + sam2 + ' > ' + pep + '_res/pred_' + pep + '_' + pdb + '_total_env1.txt'
			subprocess.call(cmd,shell=True)
		os.chdir(pep + '_res')
		envs = glob.glob('*_env1.txt')
		with open('pred_' + pep + '_total_env.txt','a') as pf :
			for env in envs :
				with open(env,'r') as ef:
					lines = ef.readlines()
					for line in lines :
						if line.find('HLA-type') < 0 :
							pf.write(line)
	#	print 'TTT'
		cmd1 = 'sort -n -k 20 pred_' + pep + '_total_env.txt > pred_' + pep + '_total_env_sort.txt' + '\n' + 'head -n 1 pred_' + pep + '_total_env_sort.txt >> ../result/pred_' + sam + '_' + sam3 + '_pre.txt'
		subprocess.call(cmd1,shell=True)
	#	print 'TTT1'
		os.chdir(wdir)

	os.chdir('result')
	# Calculation of hydrophobic count
	cmd2 = GEAR + '/hp_map pred_' + sam + '_' + sam3 + '_pre.txt > pred_' + sam + '_' + sam3 + '_pre1.txt'
	subprocess.call(cmd2,shell=True)
	
	# Annotation
	if sam3 == 'biomarker':
		cmd3 = GEAR + '/biomarker_neoscan_annot ' +  LIB + '/biomarker_5_hlatype_ntis1.txt pred_' + sam + '_' + sam3 + '_pre1.txt > pred_' + sam + '_' + sam3 + '.txt'
	elif sam3 == 'target':
		cmd3 = GEAR + '/target_neoscan_annot ' +  LIB + '/total_strong_binder.txt pred_' + sam + '_' + sam3 + '_pre1.txt > pred_' + sam + '_' + sam3 + '.txt'
	elif sam3 == 'TP53':
		cmd3 = GEAR + '/neoscan_annot ' + LIB + '/TP53_TCGA_mutation_tiled_pep.txt pred_' + sam + '_' + sam3 + '_pre1.txt > pred_' + sam + '_' + sam3 + '.txt'		
	elif sam3 == 'KRAS':
		cmd3 = GEAR + '/neoscan_annot ' + LIB + '/KRAS_TCGA_mutation_tiled_pep.txt pred_' + sam + '_' + sam3 + '_pre1.txt > pred_' + sam + '_' + sam3 + '.txt'
	subprocess.call(cmd3,shell=True)
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
