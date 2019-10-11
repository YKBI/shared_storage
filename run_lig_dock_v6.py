#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import string
import pandas as pd
import random
import glob

def pre_inp(sam,sam1):
	try:
		if not os.path.exists(sam):
			os.mkdir(sam,0777)
	except OSError:
		pass
	with open(sam + '/' + sam + '_protein.pdb','a') as pf :
		with open(sam + '.pdb','r') as tf :
			lines = tf.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 and line.find('OXT') < 0 :
					pf.write(line)
		pf.write('TER\n')
	with open(sam + '/' + sam + '_ligand.pdb','a') as lf :
		with open(sam + '.pdb','r') as tf :
			lines = tf.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 and line.find(sam1) > 0 and line.find('OXT') < 0:
					lf.write(line)
		lf.write('END\n')
	
	cmd = GEAR + '/reduce -Trim ' + sam + '/' + sam + '_ligand.pdb > ' + sam + '/' + sam + '_ligand_red.pdb'
	subprocess.call(cmd,shell=True)
	cmd1 = 'babel -ipdb ' + sam + '/' + sam + '_ligand_red.pdb -osdf ' + sam + '/' + sam + '_ligand.sdf'
	subprocess.call(cmd1,shell=True)
	cmd2 = 'babel -ipdb ' + sam + '/' + sam + '_ligand_red.pdb -omol2 ' + sam + '/' + sam + '_ligand_red.mol2'
	subprocess.call(cmd2,shell=True)


def gen_comp(sam,sam1,opt) :
	if not os.path.exists(sam + '_receptor.pdb'):
		cmdx1 = 'reduce -Trim ' + sam + '_protein.pdb > ' + sam + '_receptor.pdb'
		subprocess.call(cmdx1,shell=True)
	if opt == 'nat':
		with open('crystal_complex.pdb','a') as npdb :
			with open(sam + '_receptor.pdb','r') as rpdb :
				lines = rpdb.readlines()
				for line in lines:
					if line.startswith('ATOM') > 0 and line.find('OXT') < 0 and line.find('HOH') < 0 :
						npdb.write(line)
			npdb.write('TER\n')
			with open('%s_org.pdb'%(sam1),'r') as lpdb :
				lines = lpdb.readlines()
				for line in lines:
					if line.startswith('ATOM') > 0 :
						if line.replace('ATOM  ','HETATM').find('OXT') < 0 :
							TEXT = line.replace('ATOM  ','HETATM')[:21] + 'X' + line.replace('ATOM  ','HETATM')[22:]
							npdb.write(TEXT)
					elif line.startswith('HETATM') > 0 and line.find('OXT') < 0:
						TEXT = line[:21] + 'X' + line[22:]
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

def gen_conf(sam,sam1):
	cmd3 = ROSETTA_PY + '/molfile_to_params.py -n ' + sam1 + ' -p ' + sam1 + ' --no-pdb --keep-names ' + sam + '_ligand_red.mol2'
	subprocess.call(cmd3,shell=True)
	if not os.path.exists(sam1 + '.params'):
		cmd3 = ROSETTA_PY + '/molfile_to_params.py -n ' + sam1 + ' -p ' + sam1 + ' --no-pdb --keep-names --recharge=0 ' + sam + '_ligand_red.mol2'	
		subprocess.call(cmd3,shell=True)
	shutil.copy(sam + '_ligand_red.pdb', sam1+ '_org.pdb')
	with open(sam1 + '.params','a') as pf:
		pf.write('PDB_ROTAMERS %s_conformers.pdb\n'%(sam1))
	cmd4 = GEAR + '/balloon -f ' + LIB + '/MMFF94.mff --nconfs 50 --nGenerations 300 --rebuildGeometry ' + sam + '_ligand.sdf conformer.sdf'
	subprocess.call(cmd4,shell=True)
	cmd5 = ROSETTA_PY + '/molfile_to_params.py -n ' + sam1 + ' -p ' + sam1 + ' --no-param conformer.sdf'
	subprocess.call(cmd5,shell=True)
	atm = []
	with open(sam + '_ligand_red.pdb','r') as sf:
		lines = sf.readlines()
		for line in lines:
			if line.startswith('HETATM') > 0 :
				atm.append(line[12:16])

	pdbs = glob.glob(sam1 +'_*.pdb')
	for i in range(1,len(pdbs)):
		idx = 0
		with open(sam1 + '_' + str('%4d'%(i)).replace(' ','0') + '_c.pdb','a') as ccf:	
			with open(sam1 + '_' + str('%4d'%(i)).replace(' ','0') + '.pdb','r') as crf:
				lines = crf.readlines()
				for line in lines:
					if line.startswith('HETATM') > 0 and line[77]!='H':
						TEXT = line[:12] + atm[idx] + line[16:]	
						idx = idx + 1
						ccf.write(TEXT)
	confs = []
	confs = glob.glob(sam1 +'_*_c.pdb')
	with open(sam1 + '_conformers.pdb','a') as cf :
		for conf in confs:
			with open(conf,'r') as fx :
				lines = fx.readlines()
				for line in lines :
					cf.write(line)

def bcl_gen_conf(sam):
	cmd = BCL_GEAR + '/bcl.exe molecule:ConformerGenerator -ensemble_filenames ' + sam + '_ligand.sdf -conformers_single_file conformer.sdf'
	subprocess.call(cmd,shell=True)
	cmd1 = ROSETTA_PY + '/molfile_to_params.py -n KTL -p KTL --conformers-in-one-file conformer.sdf'
	subprocess.call(cmd1,shell=True)

sam = sys.argv[1] # PDB
sam1 = sys.argv[2] # Lig
sam2 = sys.argv[3] # No.of conformer
sam3 = sys.argv[4] # N.of CPU for rosetta run
out_arg = sys.argv[5] #log's location
DOCK_RES = 'dock_res'
ROSETTA_BIN = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/source/bin'
ROSETTA_DB = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/database'
ROSETTA_PY = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/source/scripts/python/public'
BCL_GEAR = '/awork06-1/bcl-3.6.1-Linux-x86_64'
GEAR = '/awork06-1/neoscan_gear'
LIB = '/awork06-1/neoscan_lib'

wdir = os.getcwd()

pre_inp(sam,sam1)

os.chdir(sam)

gen_conf(sam,sam1)
gen_comp(sam,sam1,'nat')
try:
	if not os.path.exists(DOCK_RES):
		os.mkdir(DOCK_RES ,0777)
except OSError:
	pass
try :
	if not os.path.exists('/awork05-1/YKLee/' + out_arg):
		os.makedirs('/awork05-1/YKLee/' + out_arg)
except OSError: pass
if os.path.exists('options'):
	os.remove('options')

with open('options','a') as of :
	of.write('-in:file:s crystal_complex.pdb\n')	
	of.write('-in:file:extra_res_fa %s.params\n'%(sam1))
	of.write('-out:path:all ' + DOCK_RES + '\n')
	of.write('-out:file:scorefile score_' + sam + '.sc\n')
	of.write('-nstruct ' + sam2 + '\n')
	of.write('-packing:ex1\n')
	of.write('-packing:ex2\n')
	of.write('-packing:no_optH false\n')
	of.write('-packing:flip_HNQ true\n')
	of.write('-packing:ignore_ligand_chi true\n')
	of.write('-parser:protocol ' + LIB + '/dock.xml\n')
	of.write('-mistakes:restore_pre_talaris_2013_behavior true\n')
	of.write('-analytic_etable_evaluation true')

if sam3 == 'single' :
	cmdx = ROSETTA_BIN + '/rosetta_scripts.linuxgccrelease @options > /awork05-1/YKLee/' + out_arg + '/' + '_'.join(sam.split('_')[:3]) + '-' + sam1 + '.run_log'
else:
	cmdx = 'mpirun -np ' + sam3 + ' ' + ROSETTA_BIN + '/rosetta_scripts.mpi.linuxgccrelease @options > /awork05-1/YKLee/' + out_arg + '/' + '_'.join(sam.split('_')[:3]) + '-' + sam1 + '.run_log'

subprocess.call(cmdx,shell=True)

with open('total_score.tsv','a') as fx :
	with open(DOCK_RES + '/score_' + sam + '.sc','r') as f :
		lines = f.readlines()
		for line in lines:
			if line.startswith('SCORE:') > 0 :
				line = "\t".join(line.split())
				fx.write(line)
				fx.write('\n')

df = pd.read_csv('total_score.tsv',sep='\t')
df[['description','total_score','total_score_X','interface_delta_X','ligand_rms_no_super_X','ligand_rms_with_super_X']].to_csv('total_score_r.tsv',sep='\t',index=False)
os.chdir('../')
