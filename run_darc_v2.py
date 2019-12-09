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

def prep_rec(sam):
	shutil.copy(wdir + '/' + sam + '.pdb',sam + '.pdb')
	with open(sam + '_rec.pdb','w') as rsf :
		with open(sam + '.pdb','r') as ssf:
			lines = ssf.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
					rsf.write(line)
	cmd = ROSETTA_BIN + '/score.default.linuxgccrelease -in:file:s ' + sam + '_rec.pdb -out:output -no_optH false'
#	print cmd
	subprocess.call(cmd,shell=True)

def make_ray(sam):
	with open(sam + '.pdb','r') as f:
		lines = f.readlines()
		for line in lines:
			if line.startswith('ATOM') > 0 :
				rch = line[21]
			elif line.startswith('HETATM') > 0 :
				olig = line[17:20]
	cmd = GEAR + '/enva_rec3 -e ' + sam + '.pdb ' + rch
	subprocess.call(cmd,shell=True)
	sips = []
	ips = []
	rips = []
	idist = []
	rips = []
	ip_df = pd.DataFrame()
	ip_df1 = pd.DataFrame()
	envs = glob.glob('*.pdb.env')
	with open(envs[0],'r') as ef:
		lines = ef.readlines()
		for line in lines:
			if line.startswith('ATOM') > 0 :
				envs = ' '.join(line[56:].split()).split(' ')
				if olig == envs[20]:
					ips.append(line[22:26].strip())
					idist.append(float(envs[19]))
		#		if float(envs[19]) < 3.5 and float(envs[19]) > 0 and olig == envs[20]:
		#			ips.append(line[22:26].strip())			
	ip_df['res'] = ips 
	ip_df['dist'] = idist
	ip_df1 = ip_df.sort_values(['dist'], ascending=[True])
	ip_df1.to_csv('ips.txt',sep='\t',index=False)
#	print ip_df1
	rips = ip_df1['res'].tolist()
#	print rips 
	sips.append(rips[0])
	sips.append(rips[1])

	with open(sam + '_het.pdb','w') as f1 :
		 with open(sam + '.pdb','r') as f2:
			lines = f2.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 :
					f1.write(line)

	cmd2 = 'babel -ipdb ' + sam + '_het.pdb -omol2 ' + sam + '_het.mol2'
	subprocess.call(cmd2,shell=True)
	lig = olig#sam.split('_')[1]
	cmd3 = ROSETTA_PY + '/molfile_to_params.py -n ' + lig + ' -p ' + lig + ' --no-pdb --keep-names ' + sam + '_het.mol2'
	subprocess.call(cmd3,shell=True)
	cmd4 = ROSETTA_BIN + '/make_ray_files.default.linuxgccrelease -pocket_static_grid -protein ' + sam + '_rec_0001.pdb -central_relax_pdb_num ' + ','.join(sips) + ' -darc_shape_only -bound_ligand ' + sam + '_het.pdb -extra_res_fa ' + lig + '.params -lig_grid -pocket_static_grid true -round_pocketGrid_center false -multiple_origin' 
#	print cmd4
	subprocess.call(cmd4,shell=True)
	return sips

def prep_lig(sam1):
	shutil.copy(wdir + '/' + sam1 + '.pdb',sam1 + '.pdb')
	cmd = GEAR + '/reduce -Trim ' + sam1 + '.pdb > ' + sam1 + '_red.pdb'
	subprocess.call(cmd,shell=True)	
	cmd2 = 'babel -ipdb ' + sam1 + '_red.pdb -omol2 ' + sam1 + '_red.mol2'
	subprocess.call(cmd2,shell=True)
	cmd4 = GEAR + '/balloon -f ' + LIB + '/MMFF94.mff --nconfs 100 --nGenerations 300 --rebuildGeometry ' + sam1 + '_red.mol2 conformer.mol2 --noGA'
	subprocess.call(cmd4,shell=True)
	with open('molfile_list.txt','w') as mf:
		mf.write('conformer.mol2\n')
	cmd3 = ROSETTA_PY + '/batch_molfile_to_params.py -d ' + ROSETTA_DB + ' --script_path=' + ROSETTA_PY + '/molfile_to_params.py molfile_list.txt' 
#	print cmd3
	subprocess.call(cmd3,shell=True)
		
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

if len(sys.argv) == 1:
	help()
else:
	sam = sys.argv[1] # PDB
	sam1 = sys.argv[2] # Lig
	sam2 = sys.argv[3] # cpu or gpu or mpi
	if sam2 == 'mpi' or sam2 =='gpu' :
		sam3 = sys.argv[4] # N.of CPU for rosetta run or gpu id

	DOCK_RES = 'dock_res'
	# Updated rosetta
	ROSETTA_BIN = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/source/bin'
	ROSETTA_DB = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/database'
	ROSETTA_PY = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/source/scripts/python/public'
	# Old rosetta 
	# ROSETTA_BIN = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/source/bin'
	# ROSETTA_DB = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/database'
	# ROSETTA_PY = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/source/scripts/python/public'
	BCL_GEAR = '/awork06-1/bcl-3.6.1-Linux-x86_64'
	GEAR = '/awork06-1/neoscan_gear'
	LIB = '/awork06-1/neoscan_lib'

	wdir = os.getcwd()

	try:
		if not os.path.exists(sam + '-' + sam1):
			os.mkdir(sam + '-' + sam1,0777)
	except OSError:
		pass

	os.chdir(sam + '-' + sam1)
	prep_rec(sam)
	ips = make_ray(sam)
	prep_lig(sam1)

	with open('options','a') as of :
		of.write('-protein %s_rec_0001.pdb\n'%(sam))
		of.write('-ligand params/conformer/000_conformers.pdb\n')
		of.write('-extra_res_fa params/conformer/000.params\n')
		of.write('-ray_file ray_%s_rec_0001_%s.txt\n'%(sam,','.join(ips)))
		of.write('-espGrid_file\n')
		of.write('-darc_shape_only\n')
		of.write('-num_particles 500\n')
		of.write('-num_runs 500\n')
		of.write('-search_conformers true\n')
		of.write('-minimize_output_complex\n')
		if sam2 == 'gpu' :
			of.write('-gpu 1\n')

#	cmdx = 'mpirun -np ' + sam2 + ' ' + ROSETTA_BIN + '/DARC.mpi.linuxgccrelease @options > run_log'
	if sam2 == 'gpu' :
		os.environ['CUDA_VISIBLE_DEVICES'] = sam3
	#	cmdz = '''export CUDA_VISIBLE_DEVICES='%d\''''%(int(sam3)) 
	#	print cmdz
	#	subprocess.call(cmdz,shell=True)
		cmdx = ROSETTA_BIN + '/DARC.opencl.linuxgccrelease @options > run_log'
	elif sam2 == 'mpi' :
		cmdx = 'mpirun -np ' + sam3 + ' ' + ROSETTA_BIN + '/DARC.mpi.linuxgccrelease @options > run_log'	
	else:
		cmdx = ROSETTA_BIN + '/DARC.linuxgccrelease @options > run_log'
	subprocess.call(cmdx,shell=True)
