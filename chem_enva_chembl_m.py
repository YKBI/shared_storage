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

def acc_chem_count(sam,lfeats,rfeats,ref):

	if not os.path.exists('crystal_complex_aa.out') :
		cmdz = GEAR + '/enva_rec3 -a crystal_complex.pdb > crystal_complex_aa.out'
		subprocess.call(cmdz,shell=True)

	mt_raccs = []
	for i in range(len(rfeats)):
		with open('crystal_complex_aa.out','r') as cc:
			lines = cc.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
					aenvs = ' '.join(line[26:].split()).split(' ')
					if rfeats[i].split('_')[0] == line[22:26].strip() and rfeats[i].split('_')[1] == line[17:20].strip() and rfeats[i].split('_')[2] == line[12:16].strip() :
						mt_raccs.append(aenvs[5])			

	tot_mt_racc = 0
	a_mt_racc = 0
	for mt_racc in mt_raccs:
		a_mt_racc = float(mt_racc) + 1
		tot_mt_racc = tot_mt_racc + a_mt_racc						

	os.chdir('dock_res')

	if not os.path.exists('ligand_feat.txt'):
		with open('ligand_feat.txt','w') as lf:
			for lfeat in lfeats:	
				if lfeat != '-' :
					lf.write('%s\n'%(lfeat))

	if not os.path.exists('feat_rmsd.txt'):
		with open('feat_rmsd.txt','a') as ff:
			ff.write('PDB\tfrmsd\tfrmsd^2\n')	
	
	idx1 = 0
	pdbs = sorted(glob.glob('*.pdb'))
	atoms = []
	for pdb in pdbs:
		print pdb
		acc = []
		phi = []
		psi = []
		sig = []
		if not os.path.exists(pdb.split('.')[0] + '.out'):
			cmd = GEAR + '/enva_rec3 -c ' + pdb.split('.')[0] + '.pdb > ' + pdb.split('.')[0] + '.out'
			subprocess.call(cmd,shell=True)
		if not os.path.exists(pdb.split('.')[0] + '_bb.out'):
			cmd1 = GEAR + '/enva_rec3 -b ' + pdb.split('.')[0] + '.pdb > ' + pdb.split('.')[0] + '_bb.out'
			subprocess.call(cmd1,shell=True)
		if not os.path.exists(pdb.split('.')[0] + '_aa.out'):
			cmd3 = GEAR + '/enva_rec3 -a ' + pdb.split('.')[0] + '.pdb > ' + pdb.split('.')[0] + '_aa.out'
			subprocess.call(cmd3,shell=True)
		if not os.path.exists('feat_rmsd.txt'):
			cmd2 = GEAR + '/rmsd_het ligand_feat.txt ' + ref + ' ' + pdb.split('.')[0] + ' >> feat_rmsd.txt'
			subprocess.call(cmd2,shell=True)
		
		raccs = []
		for i in range(len(rfeats)):
			with open(pdb.split('.')[0] + '_aa.out','r') as af :
				lines = af.readlines()
				for line in lines:
					if line.startswith('ATOM') >0 :
						aenvs = ' '.join(line[26:].split()).split(' ')
						if rfeats[i].split('_')[0] == line[22:26].strip() and rfeats[i].split('_')[1] == line[17:20].strip() and rfeats[i].split('_')[2] == line[12:16].strip() :
							raccs.append(aenvs[5])

		tot_racc = 0
		a_tot_racc = 0 
		for racc in raccs:
			a_tot_racc = float(racc) + 1
			tot_racc = tot_racc + a_tot_racc

		del_racc = 0
		del_racc = tot_racc - tot_mt_racc

		tidx = ' '
		tot_lig_acc = 0
		ave_lig_acc = 0
		tot_lig_num = 0
		num = 0
		ratio = 0
		sidx = []
		zdx = 0
		with open(pdb.split('.')[0] + '.out','r') as ef :
			lines = ef.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 :
					zdx = 0
					envs = ' '.join(line[56:].split()).split(' ')
					if int(envs[8]) == 1 :
						tot_lig_acc = tot_lig_acc + float(envs[0])
						tot_lig_num = tot_lig_num + 1
						tidx = envs[3] + '_' + envs[5] + '_' + envs[6]
						if tidx in rfeats :
							zdx = rfeats.index(tidx) + 1
							sidx.append(zdx)
							num = num + 1
		sidx.sort()
		ratio = float(num)/float(len(rfeats))
		if tot_lig_num > 0 :
			ave_lig_acc = tot_lig_acc/float(tot_lig_num)
		else : 
			ave_lig_acc = 0 

		if not os.path.exists(pdb.split('.')[0] + '.ser'):
			with open(pdb.split('.')[0] + '.ser','w') as ff4:
				for sid in sidx:
					ff4.write('%s\n'%(sid))

		cmd3 = GEAR + '/iskew ' + pdb.split('.')[0] + '.ser >> ../' + sam + '_energy_matrix/total_sk.txt'
		subprocess.call(cmd3,shell=True)
		
		with open('../' + sam + '_energy_matrix/total_sk1.txt','w') as ff5:
			if idx1 == 0 :
				ff5.write('PDB\tskewness\tClass\tDecision\n')
			with open('../' + sam + '_energy_matrix/total_sk.txt','r') as ff6:
				lines = ff6.readlines()
				for line in lines:
					if line.split('\t')[1] == 'nan':
						ff5.write(line.replace('nan','1.000000'))
					else:
						ff5.write(line)

		with open('../' + sam + '_energy_matrix/total_rac_ct.txt','a') as rf :
			if idx1 == 0 :
				rf.write('PDB')
				if len(rfeats_atm) > 0 :
					for rfeat_atm in rfeats_atm:
						rf.write('\tAA_%s'%(rfeat_atm))
				else:
					for rfeat in rfeats:
						rf.write('\tAA_%s'%(rfeat))
				rf.write('\ttotal_rec_acc\tdelta_racc\tave_lig_acc\t%Match\n')
			for i in range(len(rfeats)):
				if i == 0 :
					rf.write('%s\t%f'%(pdb.split('.')[0],float(raccs[i])+1))	
				else:
					rf.write('\t%f'%(float(raccs[i])+1))
			rf.write('\t%f\t%f\t%f\t%f\n'%(tot_racc,del_racc,ave_lig_acc,ratio))

		n_bb = 0 
		with open('../' + sam + '_energy_matrix/total_hh_ct.txt','a') as hf :
			if idx1 == 0 :
				hf.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
			with open(pdb.split('.')[0] + '_bb.out','r') as hh :
				hlines = hh.readlines()
				for hline in hlines:
					for i in range(len(rfeats)):
						if rfeats[i].split('_')[0] == hline[22:26].strip() and rfeats[i].split('_')[1] == hline[17:20].strip() and rfeats[i].split('_')[2] == hline[12:16].strip() :
							n_bb = n_bb + 1
			hf.write('%s\t%d\t%d\n'%(pdb.split('.')[0],len(hlines),n_bb))
		idx1 = idx1 + 1

	with open('../' + sam + '_energy_matrix/total_sk1.txt','w') as ff5:
		ff5.write('PDB\tskewness\tClass\tDecision\n')
		with open('../' + sam + '_energy_matrix/total_sk.txt','r') as ff6:
			lines = ff6.readlines()
			for line in lines:
				if line.split('\t')[1] == 'nan':
					ff5.write(line.replace('nan','1.000000'))
				else:
					ff5.write(line)					

	os.chdir('../')

def rmsd_calc(ref,lfeats):
	os.chdir('dock_res')
	
	pdbs = glob.glob('*_red.pdb')	

	with open('ligand_feat.txt','w') as lf:
		for lfeat in lfeats:
			if lfeat != '-' :
				lf.write('%s\n'%(lfeat))

	with open('feat_rmsd.txt','a') as ff:
		ff.write('PDB\trmsd\n')
	for pdb in pdbs:
		cmd = GEAR + '/rmsd_het ligand_feat.txt ' + ref + ' ' + pdb.split('.')[0] + ' >> feat_rmsd.txt'
		subprocess.call(cmd,shell=True)

	with open('feat_rmsd1.txt','a') as ff1:
		with open('feat_rmsd.txt','r') as ff:
			lines = ff.readlines()
			for line in lines:
				if line.startswith('PDB') > 0 :
					ff1.write('PDB\trmsd\tfrmsd\n')
				else:
					ff1.write('%s\t%s\t%s'%(line.split('\t')[0].replace('_red',''),line.split('\t')[1],line.split('\t')[2]))
	os.chdir('../')
def help():
        print "print help usage\n"
        print "Usage: python %s -s [ sample name ] -i [ input folder] -r [ ref version ] -f [ format ] -d [ cancer or non-cancer]\n"%sys.argv[0]
        return

if len(sys.argv) == 1 :
	help()
else :
	sam = sys.argv[1]
	sam1 = sys.argv[2]
	feat = ['Energy_all[bond]','Energy_all[angle]','Energy_all[dih]','Energy_all[total]','Energy_all[vdw]','Energy_all[elec]','Energy_all[elec14]','Energy_rec_lig[elec14]','Energy_rec_lig[elec]']
	header = ['PDB','Energy_all[bond]','Energy_all[angle]','Energy_all[dih]','Energy_all[vdw14]','Energy_all[elec14]','Energy_all[vdw]','Energy_all[elec]','Energy_all[total]','Energy_rec_lig[bond]','Energy_rec_lig[angle]','Energy_rec_lig[dih]','Energy_rec_lig[vdw14]','Energy_rec_lig[elec14]','Energy_rec_lig[vdw]','Energy_rec_lig[elec]','Energy_rec_lig[total]','LIE_all[EELEC]','LIE_all[EVDW]','LIE_charge[EELEC]','LIE_charge[EVDW]','LIE_nonpolar[EELEC]','LIE_nonpolar[EVDW]','LIE_polar[EELEC]','LIE_polar[EVDW]']
	cl = ['energy_all','energy_rec_lig','lie_all','lie_charge','lie_nonpolar','lie_polar']
	GEAR = '/home/user1/neoscan_gear' #'awork06-1/neoscan_gear'

	if len(sam.split('_')) > 2 :
		lig = sam.split('_')[1].split('-')[0] #sam.split('_')[len(sam.split('_'))-1]
		rec = sam.split('_')[0] #sam.split('_')[0] + '_' + sam.split('_')[1]
	else :
		lig = sam.split('_')[1].split('-')[0] #sam.split('_')[1]
		rec = sam.split('_')[0] #sam.split('_')[0]

	lfeats=[]
	rfeats=[]
	rfeats_atm=[]
	rfeats_atmn=[]
	rfeats_res=[]
	rfeats_resn=[]
	df = pd.read_csv(sam1,sep='\t')
	col_lists = df.columns.tolist()
	lfeats = df[lig].tolist()
	if 'index' in col_lists:
		rfeats = df['index'].tolist()
	else:
		rfeats_atm = df['rec_atom'].tolist()
		rfeats_atmn = df['rec_atomn'].tolist()
		rfeats_res = df['rec_res'].tolist()
		rfeats_resn = df['rec_resn'].tolist()
		for i in range(len(rfeats_res)):
			rfeats.append(str(rfeats_res[i]) + '_' + rfeats_resn[i] + '_' + rfeats_atmn[i])

#	print rfeats

	os.chdir(sam)

	try: 
		if not os.path.exists(sam + '_energy_matrix'):
			os.mkdir(sam + '_energy_matrix',0777)
	except OSError:
		pass

	acc_chem_count(sam,lfeats,rfeats,'../crystal_complex')

	with open('total_score_r1.tsv','w') as tf :
		with open('total_score_r.tsv','r') as fx:
			lines = fx.readlines()
			for line in lines:
				if line.startswith('description') > 0 :
					tf.write(line.replace('description','PDB'))
				else:
					tf.write(line)

	df_or = pd.read_csv('total_score_r1.tsv',sep='\t')
	df_sk = pd.read_csv(sam + '_energy_matrix/total_sk1.txt',sep='\t')
	df_rac = pd.read_csv(sam + '_energy_matrix/total_rac_ct.txt',sep='\t')
	df_hh = pd.read_csv(sam + '_energy_matrix/total_hh_ct.txt',sep='\t')
	df_rmsd = pd.read_csv('dock_res/feat_rmsd1.txt', sep='\t')
	total_df1 = [df_or,df_rmsd,df_hh,df_rac,df_sk]
	df_final1 = reduce(lambda left,right: pd.merge(left,right, on=['PDB'], how='outer'), total_df1)
	df_final1.to_csv(sam + '_energy_matrix/full_rec_env.txt',sep='\t',index=False,na_rep='-')
#	df_sc = pd.read_csv('total_score_r.tsv',sep='\t')
#	df_total_env = pd.merge(df_sc,df_hh)
#	df_total_env1 = pd.merge(df_total_env,df_ac)
#	df_total_env1.to_csv(sam + '_energy_matrix/full_env.txt',sep='\t',index=False,na_rep='-')
	os.chdir('../')
