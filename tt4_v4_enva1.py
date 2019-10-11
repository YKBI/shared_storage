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

def hbond_count(tag,aaa,idx) :
	flists=[]
	if os.path.exists(aaa + '/' + aaa + '_pep1.pdb') :
		os.remove(aaa + '/' + aaa + '_pep1.pdb')
	if os.path.exists(aaa + '/' + aaa + '_rec3.pdb') :
		os.remove(aaa + '/' + aaa + '_rec3.pdb')
	with open(aaa + '/' + aaa + '_pep1.pdb','a') as f1 :
		with open(aaa + '/' + aaa + '_pep.pdb','r') as f2 :
			lines1 = f2.readlines()
			for line1 in lines1 :
				if line1.startswith('ATOM') > 0 :
					text = line1.replace('ATOM  ','HETATM')
					f1.write(text)
				else :
					f1.write(line1)
	flists.append(aaa + '/' + aaa + '_rec1.pdb')
	flists.append(aaa + '/' + aaa + '_pep1.pdb')
	with open(aaa + '/' + aaa + '_rec3.pdb','a') as f3:
		for flist in flists :
			with open(flist,'r') as inf :
				f3.write(inf.read())
	cmd = GEAR + '/enva -b ' + aaa + '/' + aaa + '_rec3.pdb ' + '> ' + tag + '_energy_matrix/' + aaa + '_hh.txt'
	subprocess.call(cmd,shell=True)
	with open(tag + '_energy_matrix/total_hh_ct.txt','a') as hf :
		if idx == 0:
			hf.write('PDB\tN.of.BB\n')
		with open(tag + '_energy_matrix/' + aaa + '_hh.txt','r') as hh :
			hlines = hh.readlines()
		hf.write('%s\t%d\n'%(aaa,len(hlines)))
	#	os.remove('energy_matrix/' + aaa + '_hh.txt')

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

	try: 
		if not os.path.exists(sam + '_energy_matrix'):
			os.mkdir(sam + '_energy_matrix',0777)
	except OSError:
		pass

	samxs = []
	with open(sam + '.list','r') as f :
		lines = f.readlines()
		for line in lines :
			samxs.append(line[:-1])
	text = ""
	idx = 0
	for samx in samxs :
		hbond_count(sam,samx,idx)
		acc_count(sam,samx,idx)
		idx = idx + 1

	df_hh = pd.read_csv(sam + '_energy_matrix/total_hh_ct.txt',sep='\t')
	df_ac = pd.read_csv(sam + '_energy_matrix/total_ac_ct.txt',sep='\t')
	df_env = pd.merge(df_hh,df_ac)
	df_env.to_csv(sam + '_energy_matrix/total_env.txt',sep='\t',index=False)
	df_sc = pd.read_csv('score/' + sam + '_total_score.tsv',sep='\t')
	df_total_env = pd.merge(df_sc,df_env)
	df_total_env.to_csv(sam + '_energy_matrix/full_env.txt',sep='\t',index=False,na_rep='-')
