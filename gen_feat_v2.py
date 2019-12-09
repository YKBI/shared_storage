#/usr/bin/env python

import os
import sys
import getopt
import subprocess
import shutil
import string
import pandas as pd
from pandas import DataFrame
from collections import Counter
import math
#from sklearn import datasets, linear_model
#import matplotlib.pyplot as plt
#import statsmodels.api as sm
import numpy as np
import glob

def remove_duplicates(li):
	my_set = set()
	res = []
	for e in li:
		if e not in my_set:
			res.append(e)
			my_set.add(e)
	return res

def trim(tag,ene,st,ed,lenx,idx1) :
	tex = ''
#	print '%s\t%f\t%d\t%d'%(tag,ene,lenx,idx1)
	if idx1 == 1:
	#	print 'YYY'
		if ene < st :
			if  idx1 == lenx-1 :
				tex = '%s\t%f\n'%(tag,st)
			else :
				tex = '%s\t%f'%(tag,st)
		elif ene > ed :
			if  idx1 == lenx-1 :
			#	print 'CCC1'
				tex = '%s\t%f\n'%(tag,ed)
			else :
			#	print 'DDD1'
				tex = '%s\t%f'%(tag,ed)
		else :
			if idx1 == lenx-1 :
			#	print 'AAA1'
				tex = '%s\t%f\n'%(tag,ene)
			else :
			#	print 'BBB1'
				tex = '%s\t%f'%(tag,ene)
	else :
	#	print 'ZZZ'
		if ene < st :
			if idx1 == lenx-1 :
			#	print 'xxx'
				tex = '\t%f\n'%(st)
			else :
		 	#	print 'aaa'
				tex = '\t%f'%(st)
		elif ene > ed :
			if idx1 == lenx-1 :
		 	#	print 'CCC'
				tex = '\t%f\n'%(ed)
			else :
		 	#	print 'DDD'
				tex = '\t%f'%(ed)
		else :
			if idx1 == lenx-1 :
		 	#	print 'AAA'
				tex = '\t%f\n'%(ene)	
			#	print '%s'%(tex)
			else :
		 	#	print 'BBB'
				tex = '\t%f'%(ene)
	return tex

def integ_energy(tag1,tag,clx,idx):
	matrix = tag1 + "_energy_matrix/" +  clx + "_full_mat.txt"
	text = tag + "/min/ML_features/" + clx + "_1_1.dat"
	text1 = ''
	if os.path.exists(text) :
		with open(text,'r') as dat :
			lines = dat.readlines()
			for line in lines :
				if idx > 0 and line.startswith('#') :
					continue
				comps = " ".join(line.split()).split(' ')
				idx1 = 0
				for comp in comps :
					with open(matrix,'a') as mat :
						if idx == 0 and line.startswith('#') :
							if idx1 == 1  :
								if idx1 == len(comps)-1 :
									text1 = "PDB\t%s\n"%(comp)
								else :
									text1 = "PDB\t%s"%(comp)
								mat.write(text1)
							elif idx1 == len(comps)-1:
								text1 = "\t%s\n"%(comp)
								mat.write(text1)
							elif idx1 > 1 :
								text1 = "\t%s"%(comp)
								mat.write(text1)
						else :
							if clx == 'energy_all' :
								if  idx1 == 1 : # Energy_all[bond] ( 220 ~ 307.4)
									text1 = trim(tag,float(comp),220,307.4,len(comps),1)
								elif idx1 == 2 : # Energy_all[angle] ( 689.985 ~ 934.2)
									text1 = trim(tag,float(comp),689.985,934.2,len(comps),2)
								elif idx1 == 3: # Energy_all[dih] ( 1959.7 ~2120)
									text1 = trim(tag,float(comp),1959.7,2120,len(comps),3)
								elif idx1 == 5 : # Energy_all[elec14] ( 4834,6281.3)
									text1 = trim(tag,float(comp),4834,6281.3,len(comps),5)
								elif idx1 == 6 : # Energy_all[vdw] (-1540.7~-1326.4)
									text1 = trim(tag,float(comp),-1540.7,-1326.4,len(comps),6)
								elif idx1 == 7 : # Energy_all[elec] (-13345.1~-12560.1)
									text1 = trim(tag,float(comp),-13345.1,-12560.1,len(comps),7)
								elif idx1 == 8 : # Energy_all[total] (-5615.1 ~ -4091.73)
									text1 = trim(tag,float(comp),-5615.1,-4091.73,len(comps),8)
								else :
									text1 = '\t%f'%(float(comp))
								if idx1 > 0 :
									mat.write(text1)	
							elif clx == 'energy_rec_lig':
								if idx1 == 1:
									if idx1 == len(comps)-1 :
										text1 = "%s\t%f\n"%(tag,float(comp))
									else :
										text1 = "%s\t%f"%(tag,float(comp))
								elif idx1 == 5 :
									text1 = trim(tag,float(comp),328,2081.6,len(comps),5)
								elif idx1 == 7 :
									text1 = trim(tag,float(comp),-2791.9,-1275.7,len(comps),7)
								else :
									if idx1 == len(comps) -1:
										text1 = "\t%f\n"%(float(comp))
									else :
										text1 = "\t%f"%(float(comp))	
								if idx1 > 0 :
									mat.write(text1)
							else : 
								if  idx1 == 1 :
									if idx1 == len(comps)-1 :
										text1 = "%s\t%f\n"%(tag,float(comp))
									else :
										text1 = "%s\t%f"%(tag,float(comp))
								elif idx1 == len(comps)-1 :
									text1 = "\t%f\n"%(float(comp))
								elif idx1 > 1 :
									text1 = "\t%f"%(float(comp))
								if idx1 > 0 :
									mat.write(text1)
					idx1 = idx1 + 1

def hbond_count(sam,idx) :
	os.chdir('dock_res_' + str(idx))
	idx1 = 0
	pdbs = glob.glob('*_red.pdb')
	with open('../' + sam + '_energy_matrix/total_hh_ct.txt','a') as hf :
		for pdb in pdbs :
			cmd = GEAR + '/enva -b ' + pdb.split('.')[0].replace('_red','')  + '_red.pdb > ' + pdb.split('.')[0].replace('_red','') + '_' + str(idx) + '_bb.out'
			subprocess.call(cmd,shell=True)
			if idx1 == 0 and idx == 0 :
				hf.write('PDB\tN.of.BB\n')
			with open(pdb.split('.')[0].replace('_red','') + '_' + str(idx) + '_bb.out','r') as hh :	
				hlines = hh.readlines()
			hf.write('%s_%s\t%d\n'%(pdb.split('.')[0].replace('_red',''),str(idx),len(hlines)))
			idx1 = idx1 + 1
	os.chdir('../')

def acc_count(tag,aaa,idx) :
	flists =[]
	bbb = aaa.split('_')
	flists.append(aaa + '/' + aaa + '_rec1.pdb')
	flists.append(aaa + '/' + aaa + '_pep.pdb')
	if os.path.exists(aaa + '/' + aaa + '_rec4.pdb') :
		os.remove(aaa + '/' + aaa + '_rec4.pdb')
	if os.path.exists(aaa + '/' + bbb[0] + 'B' + bbb[1] + '_' + bbb[2] + '_rec4.pdb.env') :
		os.remove(aaa + '/' + bbb[0] + 'B' + bbb[1] + '_' + bbb[2] + '_rec4.pdb.env')
	with open(aaa + '/' + aaa + '_rec4.pdb','a') as f4:
		for flist in flists :
			with open(flist,'r') as inf :
				f4.write(inf.read())
	os.chdir(aaa)
	cmd = GEAR + '/enva -e ' + aaa + '_rec4.pdb B'
	subprocess.call(cmd,shell=True)
	os.chdir('../')
	acc = []
	phi = []
	psi = []
 	with open(aaa + '/' + bbb[0] + 'B' + bbb[1] + '_' + bbb[2] + '_rec4.pdb.env','r') as ef :
#	with open(aaa + '/' + bbb[0] + 'B' + bbb[1] + '_rec4.pdb.env','r') as ef :
		lines = ef.readlines()
		for line in lines :
			if line.find('chain') < 0 and line.startswith('ATOM') > 0 :
				envs = " ".join(line[56:].split()).split(' ')
				acc.append(envs[2])
				phi.append(envs[4])
				psi.append(envs[5])

	with open(tag + '_energy_matrix/total_ac_ct.txt','a') as af :
		if idx == 0 :
			af.write('PDB\tP1\tP2\tP3\tP4\tP5\tP6\tP7\tP8\tP9\tPHI1\tPHI2\tPHI3\tPHI4\tPHI5\tPHI6\tPHI7\tPHI8\tPHI9\tPSI1\tPSI2\tPSI3\tPSI4\tPSI5\tPSI6\tPSI7\tPSI8\tPSI9\n')
		af.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(aaa,acc[0],acc[1],acc[2],acc[3],acc[4],acc[5],acc[6],acc[7],acc[8],phi[0],phi[1],phi[2],phi[3],phi[4],phi[5],phi[6],phi[7],phi[8],psi[0],psi[1],psi[2],psi[3],psi[4],psi[5],psi[6],psi[7],psi[8]))

def acc_chem_count(sam,idx):
	os.chdir('dock_res_' + str(idx))
	
	idx1 = 0
	pdbs = glob.glob('*.pdb')
	atoms = []
	for pdb in pdbs:
		acc = []
		phi = []
		psi = []
		cmdx = GEAR + '/reduce -Trim ' + pdb.split('.')[0] + '.pdb > ' + pdb.split('.')[0] + '_red.pdb'
		subprocess.call(cmdx,shell=True)
		cmd = GEAR + '/enva_chem -c ' + pdb.split('.')[0] + '_red.pdb > ' + pdb.split('.')[0] + '_' + str(idx) + '.out'
		subprocess.call(cmd,shell=True)
		
		natom = 0 
		with open(pdb.split('.')[0] + '_' + str(idx) + '.out','r') as ef :
			lines = ef.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 and line[77]!='H' :
					envs = ' '.join(line[56:].split()).split(' ')
					if idx1 == 0 and idx == 0 :
						atoms.append(line[12:16].strip())
					acc.append(envs[1])
					phi.append(envs[2])
					psi.append(envs[3])
					natom = natom + 1

		with open('../' + sam + '_energy_matrix/total_ac_ct.txt','a') as af :
			if idx1 == 0 and idx == 0 :
				af.write('PDB')
				for atom in atoms:
					af.write('\tAA_%s'%(atom))
				for atom in atoms:
					af.write('\tPHI_%s'%(atom))		
				for atom in atoms:
					af.write('\tPSI_%s'%(atom))
				af.write('\n')
			for i in range(natom):
				j = i + 1
				if j == 1 :
					af.write('%s_%s\t%s'%(pdb.split('.')[0],str(idx),acc[i]))
				else :
					af.write('\t%s'%(acc[i]))
			for i in range(natom):
				j = i + 1
				af.write('\t%s'%(phi[i]))
			for i in range(natom):
				j = i + 1
				af.write('\t%s'%(psi[i]))
			af.write('\n')
		idx1 = idx1 + 1

	os.chdir('../')
					
def help():
        print "print help usage\n"
        print "Usage: python %s -s [ sample name ] -i [ input folder] -r [ ref version ] -f [ format ] -d [ cancer or non-cancer]\n"%sys.argv[0]
        return

if len(sys.argv) == 1 :
	help()
else :
	sam = sys.argv[1]
	feat = ['Energy_all[bond]','Energy_all[angle]','Energy_all[dih]','Energy_all[total]','Energy_all[vdw]','Energy_all[elec]','Energy_all[elec14]','Energy_rec_lig[elec14]','Energy_rec_lig[elec]']
	header = ['PDB','Energy_all[bond]','Energy_all[angle]','Energy_all[dih]','Energy_all[vdw14]','Energy_all[elec14]','Energy_all[vdw]','Energy_all[elec]','Energy_all[total]','Energy_rec_lig[bond]','Energy_rec_lig[angle]','Energy_rec_lig[dih]','Energy_rec_lig[vdw14]','Energy_rec_lig[elec14]','Energy_rec_lig[vdw]','Energy_rec_lig[elec]','Energy_rec_lig[total]','LIE_all[EELEC]','LIE_all[EVDW]','LIE_charge[EELEC]','LIE_charge[EVDW]','LIE_nonpolar[EELEC]','LIE_nonpolar[EVDW]','LIE_polar[EELEC]','LIE_polar[EVDW]']
	cl = ['energy_all','energy_rec_lig','lie_all','lie_charge','lie_nonpolar','lie_polar']
	GEAR = '/awork06-1/neoscan_gear'

	os.chdir(sam)

	pdbs = glob.glob('*.pdb')
	rec_feats=[]
	for pdb in pdbs:
		cmd = GEAR + '/enva_rec -c ' + pdb + ' > ' + pdb.split('.')[0] + '.out'
		subprocess.call(cmd,shell=True)
		trec_feats = []
		trec_feats_uniq = []
		with open(pdb.split('.')[0] + '.out','r') as f1:
			lines = f1.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 and line[77]!='H' :
					envs = ' '.join(line[56:].split()).split(' ')
					if int(envs[8]) > 0 :
						trec_feats.append(envs[4])
		trec_feats_uniq = remove_duplicates(trec_feats)
		rec_feats = rec_feats + trec_feats_uniq

#	print rec_feats
	cut = round(len(pdbs)*0.7)
	results = Counter(rec_feats)
	rec_aas = []
	for key in results:
		if int(results[key]) >= cut :
		#	print '%s\t%s'%(key,results[key])
			rec_aas.append(key)
	
	with open(sam + '_rec.txt','a') as r :
		r.write('rec_atom\trec_atomn\trec_res\trec_resn\tcrec_res\tCount\n')
		idx = 0 
		cnum = 0 
		ores = ""
		with open(pdbs[0],'r') as tf :
			lines = tf.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
				#	if idx == 0 :
				#		fres = line[22:26]
					if ores != line[22:26] :
						cnum = cnum + 1
					for key in results:
						if line[6:11].strip() == key and int(results[key]) >= cut :
						#	cnum = int(line[22:26].strip()) - int(fres) + 1
							r.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(key,line[12:16].strip(),line[22:26].strip(),line[17:20].strip(),str(cnum),results[key]))
					ores = line[22:26]
					idx = idx + 1
	for pdb in pdbs:
		if len(pdb.split('.')[0].split('_')) > 2 :
			lig = '_'.join(pdb.split('.')[0].split('_')[2:4])
		else:
			lig = pdb.split('.')[0].split('_')[1]
		with open(pdb.split('.')[0] + '_feat.out','a') as f3:
			f3.write('rec_atom\t%s\t%s_dist\n'%(lig,lig))
			with open(pdb.split('.')[0] + '.out','r') as f2:
				lines = f2.readlines()
				for line in lines:
					if line.startswith('HETATM') > 0 and line[77]!='H' :
						envs = ' '.join(line[56:].split()).split(' ')
						if envs[4] in rec_aas:
							f3.write('%s\t%s\t%s\n'%(envs[4],line[12:16].strip(),envs[7]))

	dfs = []
	df0 = pd.read_csv(sam + '_rec.txt',sep='\t')
	dfs.append(df0)
	idx = 1 
	for pdb in pdbs :
		dft = 'df_' + str(idx)
		dft = pd.read_csv(pdb.split('.')[0] + '_feat.out',sep='\t')
		dfs.append(dft)
		idx = idx + 1

#	print dfs

	df_final = reduce(lambda left,right: pd.merge(left,right, on=['rec_atom'], how='outer'), dfs)
	df_final.to_csv('%s_feat_atom_stat.txt'%(sam),sep='\t',index=False,na_rep='-')
	
	os.chdir('../')
