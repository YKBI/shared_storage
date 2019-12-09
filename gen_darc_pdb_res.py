#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import string
import pandas as pd
import random
import glob

GEAR = '/awork06-1/neoscan_gear'

sam = sys.argv[1]
sam1 = sys.argv[2]

rec = sam.split('_')[0] + '_' + sam.split('_')[1]
lig = sam.split('_')[2]

res = []
atm = []
hatm = []
with open(sam + '/' + rec + '.pdb','r') as ff:
	lines = ff.readlines()
	for line in lines:
		if line[12:16].strip() == 'CA' :
			res.append(line[22:26])

with open(sam + '/' + lig + '.pdb','r') as ff1 :
	lines = ff1.readlines()
	for line in lines:
		if line.startswith('HETATM') > 0 or line.startswith('ATOM') > 0 :
			hatm.append(line[12:16])
		#	print line[:-1]

atm_conv = {}
idx = 0
with open(sam + '/params/conformer/000_conformers.pdb','r') as ff2 :
	lines = ff2.readlines()
	for line in lines:
		if line.startswith('HETATM') > 0 :
			atm_conv[line[12:16]] = hatm[idx]
			idx = idx + 1
		elif line.startswith('TER') > 0 :
			break
#print res
#print atm
orn = ''
rn = ''
idx = 0 
idx1 = 0

cmd = GEAR + '/reduce -Trim ' + sam + '/mini_' + rec + '_rec_0001_000.pdb > ' + sam + '/mini_' + rec + '_rec_0001_000_red.pdb'
subprocess.call(cmd,shell=True)

with open(rec + '_' + lig + '.pdb','w') as ff2 :
	with open(sam + '/mini_' + rec + '_rec_0001_000_red.pdb','r') as ff1:
		lines = ff1.readlines()
		for line in lines:
			if line.startswith('ATOM') > 0 :
				if orn == '' :
					TEXT = line[:22] + res[idx] + line[26:]
				#	TEXT1 = TEXT[:6] + atm[idx1] + TEXT[11:]
					ff2.write(TEXT)
				elif orn != line[22:26] :
					idx = idx + 1
					TEXT = line[:22] + res[idx] + line[26:]
				#	TEXT1 = TEXT[:6] + atm[idx1] + TEXT[11:]
					ff2.write(TEXT)
				else :
				#	print line[:-1]
					TEXT = line[:22] + res[idx] + line[26:]
				#	print TEXT[:-1]
				#	TEXT1 = TEXT[:6] + atm[idx1] + TEXT[11:]
				#	print TEXT1[:-1]
					ff2.write(TEXT)
				orn = line[22:26]
			elif line.startswith('HETATM') > 0 :
			#	TEXT = line[:12] + atm_conv[line[12:16]] + line[16:]
				if len(lig) > 3 :
					TEXT1 = line[:17] + 'LIG' + line[20:]
				else:
					TEXT1 = line[:17] + lig + line[20:]
				ff2.write(TEXT1)
			elif line.startswith('TER') > 0 or line.startswith('END') > 0 :
				ff2.write(line)
			idx1 = idx1 + 1


ridx = []
tidx = ' '
df = pd.read_csv(sam1,sep='\t') 
ridx = df['index'].tolist()
cmd1 = GEAR + '/enva_rec3 -c ' + rec + '_' + lig + '.pdb > ' + rec + '_' + lig + '.out'
subprocess.call(cmd1,shell=True)
num = 0 
ratio = 0
sidx = []
zdx = 0 
with open(rec + '_' + lig + '.out','r') as ff3:
	lines = ff3.readlines()
	for line in lines:
		if line.startswith('HETATM') > 0 :
			zdx = 0
			envs = ' '.join(line[56:].split()).split(' ')
			tidx = envs[3] + '_' + envs[5] + '_' + envs[6]
			if tidx in ridx :
				zdx = ridx.index(tidx) + 1
				sidx.append(zdx)
				num = num + 1

sidx.sort()
ratio = float(num)/float(len(ridx))

with open(rec + '_' + lig + '.match','w') as ff6:
	ff6.write('PDB\t%Match\n')
	ff6.write('%s\t%s\n'%(rec + '_' + lig,str(ratio)))

with open(rec + '_' + lig + '.ser','w') as ff4:
	for sid in sidx:
		ff4.write('%s\n'%(sid))

with open(rec + '_' + lig +'_sk.out','w') as ff5:
	ff5.write('PDB\tskewness\tclass\tDecision\n')
cmd3 = GEAR + '/iskew ' + rec + '_' + lig + '.ser >> ' + rec + '_' + lig +'_sk.out' 
subprocess.call(cmd3,shell=True)

df1 = pd.read_csv(rec + '_' + lig + '.match',sep='\t')
df2 = pd.read_csv(rec + '_' + lig + '_sk.out',sep='\t')
df3 = pd.merge(df1,df2)
df3.to_csv(rec + '_' + lig + '_stat.txt',sep='\t',index=False)
os.remove(rec + '_' + lig + '.match')
os.remove(rec + '_' + lig + '.ser')
os.remove(rec + '_' + lig + '.out')
os.remove(rec + '_' + lig +'_sk.out')
