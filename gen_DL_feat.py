#/usr/bin/env python

import os
import sys
import getopt
import subprocess
import shutil
import string
import pandas as pd
from pandas import DataFrame
#from sklearn import datasets, linear_model
#import matplotlib.pyplot as plt
#import statsmodels.api as sm
import numpy as np
import glob

GEAR = '/awork06-1/neoscan_gear'

try:
	if not os.path.exists('final_res'):
		os.mkdir('final_res' ,0777)
except OSError:
	pass

sam = sys.argv[1]
sam1 = sys.argv[2]

pdbs = []
ligs = []
with open(sam,'r') as sf :
	lines = sf.readlines()
	for line in lines:
		pdbs.append(line[:-1])
		ligs.append(line[:-1].split('_')[1])
#pdbs = glob.glob('*.pdb')

df = pd.read_csv(sam1,sep='\t')
recs = []
recs = df['rec_atom'].tolist()
feats = []
feats.append('N.of.BB')
for rec in recs:
	feats.append('AA_' + str(rec))

angles = []
for rec in recs:
        angles.append('PHI_' + str(rec))
        angles.append('PSI_' + str(rec))

tot_df_list = []
idx = 0 
for pdb in pdbs:
	if len(pdb.split('_')) > 2:
		lig = '_'.join(pdb.split('_')[2:4])
	else :
		lig = pdb.split('_')[1]
	with open('final_res/' + pdb + '_full_rec_env.txt','w') as ff1:
		with open(pdb + '/' + pdb + '_energy_matrix/full_rec_env.txt','r') as ff:
			lines = ff.readlines()
			for line in lines:
				if line.startswith('PDB') > 0 :
					ff1.write(line)
				else:
					ff1.write(lig + '_' + line)
	with open('final_res/' + pdb + '_full_env.txt','w') as ff2:
		with open(pdb + '/' + pdb + '_energy_matrix/full_env.txt','r') as ff3:
			lines = ff3.readlines()
			for line in lines:
				if line.startswith('PDB') > 0 :
					ff2.write(line)
				else:
					ff2.write(lig + '_' + line)
	df = pd.read_csv('final_res/' + pdb + '_full_env.txt',sep='\t')
	df1 = df[df['ligand_rms_no_super_X']<1.0]
	df1[angles].drop_duplicates(angles).to_csv('final_res/' + pdb + '_angle',sep='\t',index=False)
	cmd = GEAR + '/ext_corr_dist_chem final_res/' + pdb + '_angle final_res/' + pdb + '_full_env.txt > final_res/' + pdb + '_full_env_dist.txt'
	subprocess.call(cmd,shell=True)
	df2 = pd.read_csv('final_res/' + pdb + '_full_env_dist.txt',sep='\t')
	df3 = df2[df2['Dist']<1.0]
	df3.to_csv('final_res/' + pdb + '_full_env_dist_1.txt',sep='\t',index=False)
#	df3['PDB'].to_csv('final_res/' + pdb + '_dist_1_list.txt',sep='\t',index=False)	
	df4 = pd.read_csv('final_res/' + pdb + '_full_rec_env.txt',sep='\t')
	df5 = pd.merge(df3['PDB'].to_frame(),df4)
	df5.to_csv('final_res/' + pdb + '_full_rec_env_dist_1.txt',sep='\t',index=False)
	with open('final_res/' + pdb + '_feat','w') as ffx:
		for feat in feats:
			ffx.write('%s\n'%(feat))
	idx = idx + 1
#tot_df = pd.concat(tot_df_list, axis=0,ignore_index=True)
#tot_df.to_csv('final_res/total_full_env.txt',sep='\t',index=False)
