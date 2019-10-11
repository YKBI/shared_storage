#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import string
import pandas as pd
import glob

def gen_native(aa,bb,cc) :
	with open(ROOT + '/' + aa + '/' + bb + '/' + cc + '/' + cc + '_native.pdb','a') as npdb :
		with open(ROOT + '/' + aa + '/' + bb + '/' + cc + '/' + cc + '_rec.pdb','r') as rpdb :
			lines = rpdb.readlines()
			for line in lines:
				npdb.write(line)
                npdb.write('TER\n')
		with open(ROOT + '/' + aa + '/' + bb + '/' + cc + '/' + cc + '_pep.pdb','r') as lpdb :
			lines = lpdb.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 and line.find('OXT') < 0:
					TEXT = line[:21] + 'B ' + line[23:]
					npdb.write(TEXT)
                npdb.write('END\n')

def sheba_run(ref,org) :
	cmd3 = GEAR + '/sheba_01 -x ' + ref + '.pdb ' + org + '.pdb'
	subprocess.call(cmd3,shell=True)
	cmd4 = GEAR + '/sheba_01 -t ' + org + '.trf ' + org + '.pdb'
	subprocess.call(cmd4,shell=True)
	shutil.move(org + '.pdb.pdb',org +'_tr.pdb')

def ext_pep(pdbx,ch) :
	with open(pdbx + '_pep.pdb','a') as f5 :
		with open(pdbx + '.pdb','r') as f6 :
			lines = f6.readlines()
			for line in lines :
				if line.startswith('ATOM') > 0 and line[21]==ch :
					f5.write(line)
		f5.write('END\n')

sam = sys.argv[1] # seq
sam1 = sys.argv[2] # pdb
sam2 = sys.argv[3] # HLA class
sam3 = sys.argv[4] # HLA type
sam4 = sys.argv[5] # N.of CPU for rosetta run
RES = sam + "_" + sam1 + "_" + sam2 + sam3
ROOT = '/awork04-2/NAVS/'
FLEXPEP_BIN = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/source/bin'
ROSETTA_DB = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/database'
GEAR = '/awork08/93_hong/NGS_ARS/Gear'

wdir = os.getcwd()

try:
	if not os.path.exists(RES):
		os.mkdir(RES,0777)
except OSError :
	pass

try:
	if not os.path.exists(RES + '/PDB_1ST'):
		os.mkdir(RES + '/PDB_1ST',0777)
except OSError :
	pass

try: 
	if not os.path.exists(RES + '/PDB_2ND'):
		os.mkdir(RES + '/PDB_2ND',0777)
except OSError :
	pass

try:
	if not os.path.exists(RES + '/DOCK_RES'):
		os.mkdir(RES + '/DOCK_RES',0777)
except OSError :
	pass

if not os.path.exists(ROOT + '/' + sam2 + '/' + sam3 + '/' + sam1 + '/' + sam1 + '_native.pdb'):
	gen_native(sam2,sam3,sam1)
#	with open(ROOT + '/' + sam2 + '/' + sam3 + '/' + sam1 + '/' + sam1 + '_native.pdb','a') as npdb :
#		with open(ROOT + '/' + sam2 + '/' + sam3 + '/' + sam1 + '/' + sam1 + '_rec.pdb','r') as rpdb :
#			lines = rpdb.readlines()
#			for line in lines:
#				npdb.write(line)
#		npdb.write('TER\n')
#		with open(ROOT + '/' + sam2 + '/' + sam3 + '/' + sam1 + '/' + sam1 + '_pep.pdb','r') as lpdb :
#			lines = lpdb.readlines()
#			for line in lines:
#				if line.startswith('ATOM') > 0 and line.find('OXT') < 0:
#					TEXT = line[:21] + 'B ' + line[23:] 
#					npdb.write(TEXT)
#		npdb.write('END\n')
		

with open(RES + "/prepack_flags",'a') as pf :
	pf.write('-s ' + sam1 + '_' + sam + '_inp.pdb\n')
	pf.write('-out:path:all ' + RES + '\n')
	pf.write('-out:prefix PPK_\n')
#	pf.write('-out:suffix test\n')
	pf.write('-out:file:scorefile score_ppk_' + sam1 + '_' + sam + '.sc\n')
	pf.write('-ex1\n')
	pf.write('-ex2aro\n')
	pf.write('-use_input_sc\n')
	pf.write('-flexpep_prepack\n')
	pf.write('-nstruct 1\n')

cmd = FLEXPEP_BIN + '/FlexPepDocking.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/prepack_flags > ' + RES + '/prepack.log'
subprocess.call(cmd,shell=True)

with open(RES + "/run_flags1",'a') as rf :
	rf.write('-s ' + RES + '/PPK_' + sam1 + '_' + sam + '_inp_0001.pdb\n')
	rf.write('-native ' + ROOT + '/' + sam2 + '/' + sam3 + '/' + sam1 + '/' + sam1 + '_native.pdb\n') 
	rf.write('-out:path:all ' + RES + '/PDB_1ST\n')
	rf.write('-out:file:scorefile score_1st_' + sam1 + '_' + sam + '.sc\n')
	rf.write('-nstruct 100\n')
	rf.write('-flexPepDocking:lowres_preoptimize\n')
	rf.write('-flexPepDocking:pep_refine\n')
	rf.write('-flexPepDocking:flexpep_score_only\n')
	rf.write('-ex1\n')
	rf.write('-ex2aro\n')

with open(RES + "/run_flags2",'a') as rf1 :
	rf1.write('-s ' + RES + '/PPK_' + sam1 + '_' + sam + '_inp_0001.pdb\n')
	rf1.write('-native ' + ROOT + '/' + sam2 + '/' + sam3 + '/' + sam1 + '/' + sam1 + '_native.pdb\n')
	rf1.write('-out:path:all ' + RES + '/PDB_2ND\n')
	rf1.write('-out:file:scorefile score_2nd_' + sam1 + '_' + sam + '.sc\n')
	rf1.write('-nstruct 100\n')
	rf1.write('-flexPepDocking:pep_refine\n')
	rf1.write('-flexPepDocking:flexpep_score_only\n')
	rf1.write('-ex1\n')
	rf1.write('-ex2aro\n')

cmd1 = 'mpirun -np ' + sam4 + ' ' + FLEXPEP_BIN + '/FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/run_flags1 > ' + RES + '/run1.log' 
cmd2 = 'mpirun -np ' + sam4 + ' ' + FLEXPEP_BIN + '/FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/run_flags2 > ' + RES + '/run2.log'
#print cmd1
#print cmd2
subprocess.call(cmd1,shell=True)
subprocess.call(cmd2,shell=True)

with open(RES + '/total_score.tsv','a') as fx :
	with open(RES + '/PDB_1ST/score_1st_' + sam1 + '_' + sam + '.sc','r') as f :
		lines = f.readlines()
		for line in lines:
			if line.startswith('SCORE:') > 0 :
				line = " ".join(line.split())
				col = line.split(' ')
				if line.find('total_score') > 0 :
					fx.write('%s\t%s\t%s\n'%(col[68],col[45],col[1]))
				else :
					fx.write('PDB_1ST/%s.pdb\t%s\t%s\n'%(col[68],col[45],col[1]))
	with open(RES + '/PDB_2ND/score_2nd_' + sam1 + '_' + sam + '.sc','r') as f1 :
		lines = f1.readlines()
		for line in lines :
			if line.startswith('SCORE:') > 0 and line.find('total_score') < 0 :
				line = " ".join(line.split())
				col = line.split(' ')
				fx.write('PDB_2ND/%s.pdb\t%s\t%s\n'%(col[62],col[45],col[1]))

df = pd.read_csv(RES + '/total_score.tsv',sep='\t')
df = df.sort_values(['total_score'])
df1 = df.head(10)
df1.to_csv(RES + '/dock_selected.tsv',sep='\t',index=False)

with open(RES + '/dock_selected.tsv','r') as f2 :
	lines = f2.readlines()
	for line in lines :
		cols = line.split('\t')
		FFS = cols[0].split('/')
		if FFS[0].find('1ST') > 0 :
			shutil.copy(RES + '/' + cols[0],RES + '/DOCK_RES/1ST_' + FFS[1])
		elif FFS[0].find('2ND') > 0 :
			shutil.copy(RES + '/' + cols[0],RES + '/DOCK_RES/2ND_' + FFS[1])

pdbs=[]

os.chdir(RES + '/DOCK_RES')
dock_list = glob.glob('*.pdb')
os.chdir(wdir)

with open(ROOT + '/' + sam2 + '_list/' + sam3 + '.list','r') as f3 :
	lines = f3.readlines()
	for line in lines :
		pdbs.append(line[:-1])

for pdb in pdbs :
	if not os.path.exists(ROOT + '/' + sam2 + '/' + sam3 + '/' + pdb + '/' + pdb + '_native.pdb'):
		gen_native(sam2,sam3,pdb)
	shutil.copy(ROOT + '/' + sam2 + '/' + sam3 + '/' + pdb + '/' + pdb + '_native.pdb',RES + '/DOCK_RES')

os.chdir(RES + '/DOCK_RES')
for pdb in pdbs :
	if pdb == sam1 :
		if not os.path.exists(pdb + '_native_pep.pdb') :
			ext_pep(pdb + '_native','B')
		for dock in dock_list :
			if not os.path.exists(dock.split('.')[0] + '_pep.pdb'):
				ext_pep(dock.split('.')[0],'B')
			cmd4 = GEAR + '/rmsd_total ' + pdb + '_native_pep.pdb ' + dock.split('.')[0] + '_pep.pdb bb >> rmsd_total.txt'
			subprocess.call(cmd4,shell=True)
	else :
		sheba_run(sam1 + '_native',pdb + '_native')
		ext_pep(pdb + '_native_tr','B')
		for dock in dock_list :
			if not os.path.exists(dock.split('.')[0] + '_pep.pdb'):
				ext_pep(dock.split('.')[0],'B')
			cmd4 = GEAR + '/rmsd_total ' + pdb + '_native_tr_pep.pdb ' + dock.split('.')[0] + '_pep.pdb bb >> rmsd_total.txt'
			subprocess.call(cmd4,shell=True)

df2 = pd.read_csv('rmsd_total.txt',sep='\t',names=('DOCK_STR','HLA-REF','PEP_RMSD'))
df2 = df2.sort_values(['PEP_RMSD'])
df2.to_csv('rmsd_total_sort.txt',sep='\t',index=False)
os.chdir(wdir)
