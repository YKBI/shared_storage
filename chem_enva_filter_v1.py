#/usr/bin/env python

import os
import sys
import getopt
import subprocess
import shutil
import string
import pandas as pd
from pandas import DataFrame
from scipy import stats
from scipy.stats import zscore
#from sklearn import datasets, linear_model
#import matplotlib.pyplot as plt
#import statsmodels.api as sm
import numpy as np
import glob
import math

GEAR = '/awork06-1/neoscan_gear'

try:
	if not os.path.exists('model'):
		os.mkdir('model' ,0777)
except OSError:
	pass

sam = sys.argv[1]

rec = sam.split('_')[0]
lig = sam.split('_')[1]

wdir = os.getcwd()

if os.path.exists(sam + '_trans/' + rec + '_' + lig):
	samx = sam + '_trans/' + rec + '_' + lig
elif os.path.exists(sam + '_trans/' + rec + '_' + lig +'_' + rec + '_' + lig + '_' + lig):
	samx = sam + '_trans/' + rec + '_' + lig +'_' + rec + '_' + lig + '_' + lig

df = pd.read_csv(samx + '/' + sam + '_energy_matrix/full_rec_env.txt',sep='\t')

cols = df.columns.tolist()

zfeats = []
for col in cols:
	if col == 'PDB' or col == '%Match' or col == 'skewness' or col == 'total_score' or col == 'total_score_X' or col == 'interface_delta_X' or col == 'ligand_rms_no_super_X' or col == 'ligand_rms_with_super_X' or col == 'N.of.BB_full' or col == 'delta_racc' or col == 'Class' or col =='Decision' or col == 'rmsd':
		continue
	zfeats.append(col)

if 'rmsd' in cols:
	df1 = df[(df['rmsd']<1.0) & (df['%Match'] >=0.4) & (df['skewness'] <=0.6)]
elif 'ligand_rms_no_super_X' in cols:
	df1 = df[(df['ligand_rms_no_super_X']<1.0) & (df['%Match'] >=0.4) & (df['skewness'] <=0.6)]
#df1 = df[df['%Match'] >=0.4]
df1.to_csv('model/' + sam + '_pre_filt.txt',sep='\t',index=False)
cond = []
with open('model/' + sam + '_cut_off.txt','a') as fz:
	fz.write('Feature\tmin\tmax\n')
for zfeat in zfeats:
	df_z = pd.DataFrame()
	df_z1 = pd.DataFrame()
	min_val = 0
	max_val = 0 
	s_zfeat_lsts = []
	zfeat_lsts = []
	f_zfeat_lsts = []
	if df1[zfeat].std() == 0 :
		if df1[zfeat].mean() == 1:
			cond.append('%s == 1'%(zfeat))
			with open('model/' + sam + '_cut_off.txt','a') as fz:
				fz.write('%s\t1\t1\n'%(zfeat))
		else:
			cond.append('%s == %d'%(zfeat,df1[zfeat].mean()))
			with open('model/' + sam + '_cut_off.txt','a') as fz:
				fz.write('%s\t%d\t%d\n'%(zfeat,df1[zfeat].mean(),df1[zfeat].mean()))
	else :
		zfeat_lsts = list(df1[zfeat])
		for zfeat_lst in zfeat_lsts:
			f_zfeat_lsts.append(math.log10(float(zfeat_lst)))
		f_zfeat_lsts = stats.zscore(f_zfeat_lsts)
		df_z[zfeat] = zfeat_lsts
		df_z[zfeat + '_zscore'] = f_zfeat_lsts
		df_z = df_z.sort_values([zfeat + '_zscore'], ascending=[True])	
	#	df_z.to_csv(zfeat + '_zs.txt',sep='\t',index=False)
		df_z1 = df_z[(df_z[zfeat + '_zscore'] > -2) & (df_z[zfeat + '_zscore'] <2)]
		s_zfeat_lsts = list(df_z1[zfeat])
		if zfeat == 'N.of.BB_feat' or zfeat == 'total_rec_acc' or zfeat == 'ave_lig_acc':
			if zfeat == 'N.of.BB_feat':
				min_bump = float(s_zfeat_lsts[0])
				max_bump = float(s_zfeat_lsts[len(s_zfeat_lsts)-1])
				with open('model/' + sam + '_cut_off.txt','a') as fz:
					fz.write('%s\t%f\t%f\n'%(zfeat,min_bump,max_bump))
			elif zfeat == 'total_rec_acc' :
				min_racc = float(s_zfeat_lsts[0])
				max_racc = float(s_zfeat_lsts[len(s_zfeat_lsts)-1])
				with open('model/' + sam + '_cut_off.txt','a') as fz:
					fz.write('%s\t%f\t%f\n'%(zfeat,min_racc,max_racc))
			elif zfeat == 'ave_lig_acc':
				min_lacc = float(s_zfeat_lsts[0])
				max_lacc = float(s_zfeat_lsts[len(s_zfeat_lsts)-1])
				with open('model/' + sam + '_cut_off.txt','a') as fz:
					fz.write('%s\t%f\t%f\n'%(zfeat,min_lacc,max_lacc))
		else:
			min_val = float(s_zfeat_lsts[0])
			max_val = float(s_zfeat_lsts[len(s_zfeat_lsts)-1])
			cond.append('%s > %f and %s < %f'%(zfeat,min_val,zfeat,max_val))
			with open('model/' + sam + '_cut_off.txt','a') as fz:
				fz.write('%s\t%f\t%f\n'%(zfeat,min_val,max_val))

#print cond
#condx = ' and '.join(cond)
#print condx

#print df
#dff = df.query(condx)
#dff.to_csv('test.txt',sep='\t',index=False)

#print '%f\t%f\t%f\t%f\t%f\t%f'%(min_bump,max_bump,min_racc,max_racc,min_lacc,max_lacc)

dff = df[(df['%Match'] >=0.4) & (df['skewness']<=0.6) & (df['N.of.BB_feat'] > min_bump)& (df['N.of.BB_feat'] < max_bump) & (df['total_rec_acc'] > min_racc) & (df['total_rec_acc'] < max_racc) & (df['ave_lig_acc'] > min_lacc) & (df['ave_lig_acc'] < max_lacc)]
dff.to_csv('model/' + sam + '_full_rec_env_temp.txt',sep='\t',index=False,na_rep='-')
#cmd = GEAR + '/acc_filter cut_off.txt temp_final.txt > final_filter.txt'
cmd = GEAR + '/acc_filter model/' + sam + '_cut_off.txt ' + 'model/' + sam + '_full_rec_env_temp.txt > model/' + sam + '_full_rec_env_f.txt'
subprocess.call(cmd,shell=True)
