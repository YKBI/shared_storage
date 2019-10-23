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
		fname = line.split('_')
		pdbs.append(line.strip())
		ligs.append(fname[-1].strip())
		#pdbs.append(line[:-1])
		#ligs.append(line[:-1].split('_')[1])
#pdbs = glob.glob('*.pdb')

tot_df_list = []
idx = 0 
for pdb in pdbs:
	lig = pdb.split('_')[-1]
	print lig
	with open('final_res/' + pdb + '_full_env.txt','w') as ff1:
		with open(pdb + '/' + pdb + '_energy_matrix/full_env.txt','r') as ff:
			lines = ff.readlines()
			for line in lines:
				if line.startswith('PDB') :
					ff1.write(line)
				else:
					ff1.write(lig + '_' + line)
	with open('final_res/total_full_env.txt','a') as ff2:
		with open('final_res/' + pdb + '_full_env.txt','r') as ff1:
			lines = ff1.readlines()
			for line in lines:
				if line.startswith('PDB') > 0 and idx == 0 :
					ff2.write(line)
				else:
					if line.startswith('PDB') > 0 :
						continue
					ff2.write(line)
		os.remove('final_res/' + pdb + '_full_env.txt')
	idx = idx + 1

cmd = GEAR + '/missing_imp final_res/total_full_env.txt > final_res/total_full_env1.txt\n sort -n -k 2 final_res/total_full_env1.txt > final_res/total_full_env1_sort.txt'
subprocess.call(cmd,shell=True)

df = pd.read_csv(sam1,sep='\t')
recs = []
recs = df['rec_atom'].tolist()

angles = []
for rec in recs:
	angles.append('PHI_' + str(rec))
	angles.append('PSI_' + str(rec))

for lig in ligs:
	with open('final_res/' + lig + '_env.txt','w') as tf :
		with open('final_res/total_full_env1_sort.txt','r') as of :
			lines = of.readlines()
			for line in lines:
				if line.startswith('PDB') > 0 or line.startswith(lig) > 0 :
					tf.write(line)
	df = pd.read_csv('final_res/' + lig + '_env.txt',sep='\t')
	df1 = df[df['rmsd']<1.0]
	df1[angles].drop_duplicates(angles).to_csv('final_res/' + lig + '_angle',sep='\t',index=False)

#tot_df = pd.concat(tot_df_list, axis=0,ignore_index=True)
#tot_df.to_csv('final_res/total_full_env.txt',sep='\t',index=False)
