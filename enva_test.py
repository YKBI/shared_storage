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
import multiprocessing
import glob

def enva_working(pdb):
	inpdb = pdb.split('.')[0]
	csplit_cmd = 'csplit -f "%s_" %s.pdb \'/TER/\';sed \'s/ATOM  /HETATM/\' %s_01 > %s_02; cat %s_00 %s_02 > %s_het.pdb'%(inpdb,inpdb,inpdb,inpdb,inpdb,inpdb,inpdb)
	subprocess.call(csplit_cmd,shell=True)
	cp_cmd = 'cp %s.pdb %s.pdb'%(inpdb,inpdb.split('_')[-3])
	subprocess.call(cp_cmd,shell=True)
	if not os.path.exists(inpdb + '_a.out'):
		os.system(GEAR + '/enva_rec.v1.1 -a ' + inpdb + '.pdb > ' + inpdb + '_a.out')
	if not os.path.exists(inpdb + '_b.out'):
		os.system(GEAR + '/enva_rec.v1.1 -b ' + inpdb + '_het.pdb > ' + inpdb + '_b.out')
	if not os.path.exists(inpdb + '_m.out'):
		os.system(GEAR + '/enva_rec.v1.1 -m ' + inpdb.split('_')[-3] + '.pdb ' + inpdb.split('_')[-3] + 'B > ' + inpdb + '_m.out')
	if not os.path.exists('*.env'):
		os.system(GEAR + '/enva_rec.v1.1 -e ' + inpdb + '.pdb B')
	cmd = 'rm -rf ' + inpdb.split('_')[-3] + '.pdb'
	subprocess.call(cmd,shell=True)
def multi_part1(x):
	pool = multiprocessing.Pool(10)
	pool.map(enva_working,x)
	pool.close()
	pool.join()
def txt_writing(x):
	ff_list = []
	hh_list = []
	tt_acc = []
	tt_phi = []
	tt_psi = []
	tt_ac_dic = {}
	mat_lst = []
	lig_acc = []
	lig_num = []
	ave_acc = []
	tot_lig_acc = 0
	ave_lig_acc = 0
	tot_lig_num = 0
	ratio = 0
	zdx = 0
	pdbl = []
	tttt = []
	newl = []
	for i in sorted(glob.glob('*rev.pdb')):
		os.system('csplit -f \'' + i + '_\' ' + i + ' \'/TER/\'')
		os.system('/awork06-1/YKLee/c_script/rmsd_total ' + PEPLIB + '/' + inpdb + '_pep.pdb ' + i + '_01 bb >> ' + x + '_rmsd.txt')
	#rdf = pd.read_csv('rmsd.txt',sep='\s+',header=None)
	for f,g in zip(sorted(glob.glob('*a.out')),sorted(glob.glob('*m.out'))):
		sidx= []
		pdb = '_'.join([f.split('.')[0].split('_')[1],str(x),f.split('.')[0].split('_')[5]])#'_'.join(f.split('.')[0].split('_')[:6])
		pdbl.append(pdb)
		adf = pd.read_csv(f,sep='\s+',skiprows=[0,0],header=None)
		mdf = pd.read_csv(g,sep='\s+',header=None)
		mdf_filter = mdf[mdf[17] == 1].iloc[:,9:16]
		for i in mdf_filter[9]:
			tot_lig_acc = tot_lig_acc + float(i)
			tot_lig_num += 1
			if tot_lig_num > 0: 
				ave_lig_acc = tot_lig_acc/float(tot_lig_num)
			else: 
				ave_lig_acc = 0
		ave_acc.append(ave_lig_acc)
		
		new_acc = {}
		num = 0
		for i in rfeats: #mdf_filter.T:
			new_acc['AA_' + cc_dic[str(i)]] =sum([adf[adf[1] == int(i.split('_')[0])][adf[2] == i.split('_')[2]][adf[3] == i.split('_')[1]][11].values,1])
		for i in mdf_filter.T:
			tidx = str(mdf_filter.T[i][13]) + '_' + str(mdf_filter.T[i][14]) + '_' + str(mdf_filter.T[i][15])
			if tidx in rfeats:
				tttt.append(tidx)
				zdx = rfeats.index(tidx) +1
				sidx.append(zdx)
				num += 1
		newl.append(pd.DataFrame(new_acc))
		if not os.path.exists(pdb.split('.')[0] + '.ser'):
			with open(pdb.split('.')[0] + '.ser','w') as W:
				for sid in list(set(sidx)):
					W.write(str(sid)+'\n')
		iskew_cmd = GEAR +'/iskew ' + pdb.split('.')[0] + '.ser >> total_' + x + '_sk.txt'
		subprocess.call(iskew_cmd,shell=True)
		ratio = float(num)/float(len(rfeats))
		mat_lst.append(ratio)
	new_ac = pd.concat(newl)
	new_ac['total_rec_acc'] = new_ac.sum(axis=1)
	new_ac['ave_lig_acc'] = ave_acc
	new_ac['%Match'] = mat_lst
	#new_ac['rmsd'] = rdf[2]
	new_ac['PDB'] = pdbl
	new_ac.set_index('PDB').reset_index()
	
	new_ac.to_csv(x + '_nac.txt',sep='\t',index=False)
	nn_list = []
	for f in sorted(glob.glob('*b.out')):
		with open(f,'r') as F:
		    n = 0
		    for line in F.readlines():
			tt = line[7:13].strip() + '_' + line[17:21].strip() + '_' + line[13:17].strip()
			if tt in rfeats:
			    n += 1
		    nn_list.append(n)
		with open(f,'r') as F:
		    hh_list.append(len(F.readlines()))
		ff_list.append(f.split('.')[0].split('_')[1]+ '_' + str(x) + '_' + f.split('.')[0].split('_')[-4])
	with open(x + '_hh.txt','w') as W:
		W.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
		for i,j,k in zip(ff_list,hh_list,nn_list):
			W.write(str(i) + '\t' + str(j) + '\t' + str(k) + '\n')		
	acc_columns = ['P%d'%i for i in range(1,seqlen+1)]
	phi_columns = ['PHI%d'%i for i in range(1,seqlen+1)]
	psi_columns = ['PSI%d'%i for i in range(1,seqlen+1)]
	pdbl = []
	for f in sorted(glob.glob('*.env')):
		pdb = '_'.join([sam.split('_')[1],str(x),f.split('.')[0].split('_')[5]])
		pdbl.append(pdb)
		with open(f,'r') as F:
			acc = []
			psi = []
			phi = []
			for line in F.readlines():
				if line.find('chain') and line.startswith('ATOM'):
					envs = line.split()[9:]
					acc.append(envs[2])
					phi.append(envs[4])
					psi.append(envs[5])
			tt_acc.append(pd.DataFrame(acc, index= acc_columns).T)
			tt_phi.append(pd.DataFrame(phi, index = phi_columns).T)
			tt_psi.append(pd.DataFrame(psi, index = psi_columns).T)
	pdbdf = pd.DataFrame(pdbl,columns = ['PDB'])
	acc_df = pd.concat(tt_acc,ignore_index=True).reindex()
	phi_df = pd.concat(tt_phi,ignore_index=True).reindex()
	psi_df = pd.concat(tt_psi,ignore_index=True).reindex()
	act_df = pd.concat([pdbdf,acc_df,phi_df,psi_df],axis=1)
	act_df.to_csv(x + '_ac_ct.txt',sep='\t',index=False)
	
#x=sam,
def work_part(x):
	os.chdir('PDB_1ST')
	os.system('rm -rf *.env *.txt *.ser *.out *rev.pdb')
	for i in sorted(glob.glob('*[0-9].pdb')):
		print i
		os.system('reduce -Trim ' + i + ' > ' + i.split('.')[0] + '_red.pdb')
		os.system('python /awork06-1/YKLee/py_script/pep_atom_revi.py ../../' + '_'.join(i.split('_')[1:4]) + '.pdb ' + i.split('.')[0] + '_red.pdb')
	print "########## End Of Revise ##########"
	lst1 = sorted(glob.glob('*rev.pdb'))
	
	multi_part1(lst1)
	cmd = 'rm -rf *_02 *_01 *_00 *_red.pdb'
	subprocess.call(cmd,shell=True)
	print "########## End Of Enva Work ##########"
	txt_writing('1ST')
	os.system('rm -rf *_00 *_01 *.ser *rev.pdb')
	ac_df1 = pd.read_csv('1ST_ac_ct.txt',sep='\t')
	nac_df1 = pd.read_csv('1ST_nac.txt',sep='\t')
	sk_df1 = pd.read_csv('total_1ST_sk.txt',sep='\t',header=None).fillna('1.000000')
	hh_df1 = pd.read_csv('1ST_hh.txt',sep='\t')
	rmsd_df1 = pd.read_csv('1ST_rmsd.txt',sep='\t',header=None)
	os.chdir('../PDB_2ND')
	os.system('rm -rf *.env *.txt *.ser *.out ')
	for i in sorted(glob.glob('*[0-9].pdb')):
                os.system('reduce -Trim ' + i + ' > ' + i.split('.')[0] + '_red.pdb')
                os.system('python /awork06-1/YKLee/py_script/pep_atom_revi.py ../../' + '_'.join(i.split('_')[1:4]) + '.pdb ' + i.split('.')[0] + '_red.pdb')
	print "########## End Of Revise ##########"
	lst2 = sorted(glob.glob('*rev.pdb'))
	multi_part1(lst2)
	subprocess.call(cmd,shell=True)
	print "########## End Of Enva Work ##########"
	txt_writing('2ND')
	os.system('rm -rf *_00 *_01 *.ser ')
	ac_df2 = pd.read_csv('2ND_ac_ct.txt',sep='\t')
	nac_df2 = pd.read_csv('2ND_nac.txt',sep='\t')
	sk_df2 = pd.read_csv('total_2ND_sk.txt',sep='\t',header=None).fillna('1.000000')
	hh_df2 = pd.read_csv('2ND_hh.txt',sep='\t')
	rmsd_df2 = pd.read_csv('2ND_rmsd.txt',sep='\t',header=None)
	os.chdir('../' + x + '_energy_matrix/')
	pd.concat([rmsd_df1,rmsd_df2]).to_csv('total_rmsd.txt',sep='\t',index=False)
	pd.concat([ac_df1,ac_df2]).to_csv('total_ac_ct.txt',sep='\t',index=False)
	pd.concat([nac_df1,nac_df2]).to_csv('total_nac.txt',sep='\t',index=False)
	pd.concat([hh_df1,hh_df2]).to_csv('total_hh_ct.txt',sep='\t',index=False)
	pd.concat([sk_df1,sk_df2]).to_csv('total_sk1.txt',sep='\t',index=False,header=['PDB','skewness','Class','Decision'])
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
	PDBLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/Native/' + sam.split('_')[2] + '_native'
	PEPLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/peptide/' + sam.split('_')[2]
	seq = sam.split('_')[0]
	inpdb = sam.split('_')[1]
	seqlen = len(seq)
	os.chdir(sam)
        try:
                if not os.path.exists(sam + '_energy_matrix'):
                        os.mkdir(sam + '_energy_matrix',0777)
        except OSError:
                pass
	rfeats = []
	feats = pd.read_csv(PDBLIB + '/' + sam.split('_')[1] + '.out',sep='\s+',header=None)
	feats_filter = feats[feats[17] ==1].iloc[:,13:16]
	for i in feats_filter.values.tolist():
		rfeats.append(str(i[0]) + '_' + '_'.join(i[1:]))
	cc_dic = {}
	feats1_filter =feats[feats[17] ==1].iloc[:,12:16]
	rfeats1 = []
	for i in feats1_filter.values.tolist():
		rfeats1.append(str(i[0]) + '_' + '_'.join(i[2:]))
	for i,j in zip(rfeats,rfeats1):
		cc_dic[i] = j
	work_part(sam)
	os.chdir('../')
	df_hh = pd.read_csv(sam + '_energy_matrix/total_hh_ct.txt',sep='\t')
	df_nac = pd.read_csv(sam + '_energy_matrix/total_nac.txt',sep='\t')
	df_sk = pd.read_csv(sam + '_energy_matrix/total_sk1.txt',sep='\t')
	df_rmsd = pd.read_csv(sam + '_energy_matrix/total_rmsd.txt',sep='\t')
	df_rmsd.columns=['pdb','reference','rmsd']
	total_df1 = [df_hh,df_nac,df_sk]
	df_final1 = reduce(lambda left,right: pd.merge(left,right, on=['PDB'],how='outer'), total_df1)
	df_final1['rmsd'] = df_rmsd['rmsd']
	df_final2 = df_final1.set_index('rmsd').reset_index().set_index('PDB').reset_index()
	df_final2.to_csv(sam + '_energy_matrix/' + sam + '_total.txt',sep='\t',index=False,na_rep='-')
	os.chdir('../')
