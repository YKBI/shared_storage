#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import string
import pandas as pd
import random
import glob
import Bio.PDB
#import PeptideBuilder
#from PeptideBuilder import Geometry
from collections import Counter

def remove_duplicates(li):
	my_set = set()
	res = []
	for e in li:
		if e not in my_set:
			res.append(e)
			my_set.add(e)
	return res

def pep_build(seq):
	if not os.path.exists('model'):
		os.mkdir('model',0777)

	with open('inp.fasta','w') as xf:
		xf.write('>input\n')
		xf.write('%s\n'%(seq))

	cmd = ROSETTA_BIN + '/BuildPeptide.default.linuxgccrelease -in:file:fasta inp.fasta -phi -139 -psi -135 -out:file:o model/peptide1.pdb'	
	subprocess.call(cmd,shell=True)
 	cmd1 = ROSETTA_BIN + '/BuildPeptide.default.linuxgccrelease -in:file:fasta inp.fasta -phi -57 -psi -47 -out:file:o model/peptide2.pdb'
#	cmd1 = ROSETTA_BIN + '/BuildPeptide.default.linuxgccrelease -in:file:fasta inp.fasta -phi -139 -psi -135 -out:file:o model/peptide2.pdb'
	subprocess.call(cmd1,shell=True)
 	cmd2 = ROSETTA_BIN + '/BuildPeptide.default.linuxgccrelease -in:file:fasta inp.fasta -phi -78 -psi 149 -out:file:o model/peptide3.pdb'
#	cmd2 = ROSETTA_BIN + '/BuildPeptide.default.linuxgccrelease -in:file:fasta inp.fasta -phi -139 -psi -135 -out:file:o model/peptide3.pdb'
	subprocess.call(cmd2,shell=True)

	with open('ens.list','w') as ef:
		ef.write('model/peptide1r.pdb\n')
		ef.write('model/peptide2r.pdb\n')
		ef.write('model/peptide3r.pdb\n')

	with open('ens-aa.list','w') as eaf:
		eaf.write('model/peptide1-aa.pdb\n')
		eaf.write('model/peptide2-aa.pdb\n')
		eaf.write('model/peptide3-aa.pdb\n')

	with open('ens-rmsd.list','w') as erf:
		erf.write('model/peptide1-heavy.pdb\n')
		erf.write('model/peptide2-heavy.pdb\n')
		erf.write('model/peptide3-heavy.pdb\n')

	

def pep_buildx(seq):
#	structure = PeptideBuilder.make_extended_structure(seq)
#	out = Bio.PDB.PDBIO()
#	out.set_structure(structure)
#	out.save('peptide.pdb')

	if not os.path.exists('model'):
		os.mkdir('model',0777)

	phi1 = []
	psi1 = []
	phi2 = []
	psi2 = []
	phi3 = []
	psi3 = []
	for i in range(len(seq)-2):
		phi1.append('-57')
		psi1.append('-47')
		phi2.append('-139')
		psi2.append('-135')
		phi3.append('-78')
		psi3.append('149')

	structure = PeptideBuilder.make_structure(seq,phi1,psi1)
	out = Bio.PDB.PDBIO()
	out.set_structure(structure)
	out.save('model/peptide1.pdb')

	structure1 = PeptideBuilder.make_structure(seq,phi2,psi2)
	out1 = Bio.PDB.PDBIO()
	out1.set_structure(structure1)
	out1.save('model/peptide2.pdb')

	structure2 = PeptideBuilder.make_structure(seq,phi3,psi3)
	out2 = Bio.PDB.PDBIO()
	out2.set_structure(structure2)
	out2.save('model/peptide3.pdb')

	with open('ens.list','w') as ef:
		ef.write('model/peptide1r.pdb\n')
		ef.write('model/peptide2r.pdb\n')
		ef.write('model/peptide3r.pdb\n')

	with open('ens-aa.list','w') as eaf:
		eaf.write('model/peptide1-aa.pdb\n')
		eaf.write('model/peptide2-aa.pdb\n')
		eaf.write('model/peptide3-aa.pdb\n')

	with epen('ens-rmsd.list','w') as erf:
		erf.write('model/peptide1-heavy.pdb\n')
		erf.write('model/peptide2-heavy.pdb\n')
		erf.write('model/peptide3-heavy.pdb\n')	
		
def pre_inp(sam,sam1):
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

def prep_rec(sam,rch,pch):
	shutil.copy(wdir + '/' + sam + '.pdb',sam + '.pdb')
	with open(sam + '_rec.pdb','w') as rsf :
		with open(sam + '.pdb','r') as ssf:
			lines = ssf.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 and line[21] == rch:
					rsf.write(line)
		rsf.write('END\n')
	
	with open(sam + '_pep.pdb','w') as csf:
		with open(sam + '.pdb','r') as ssf:
			lines = ssf.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 and line[21] == pch:
					csf.write(line)
		csf.write('END\n')

def make_ray(sam):
	with open(sam + '_conv.pdb','r') as f:
		lines = f.readlines()
		for line in lines:
			if line.startswith('ATOM') > 0 :
				rch = line[21]
			elif line.startswith('HETATM') > 0 :
				olig = line[17:20]
	cmd = GEAR + '/enva_rec3 -e ' + sam + '_conv.pdb ' + rch
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
		with open(sam + '_conv.pdb','r') as f2:
			lines = f2.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 :
					f1.write(line)

	cmd2 = 'babel -ipdb ' + sam + '_het.pdb -omol2 ' + sam + '_het.mol2'
	subprocess.call(cmd2,shell=True)
	cmd3 = ROSETTA_PY + '/molfile_to_params.py -n pep -p pep --no-pdb --keep-names ' + sam + '_het.mol2'
	subprocess.call(cmd3,shell=True)
	cmd4 = ROSETTA_BIN + '/make_ray_files.default.linuxgccrelease -pocket_static_grid -protein ' + sam + '_rec_0001.pdb -central_relax_pdb_num ' + ','.join(sips) + ' -darc_shape_only -bound_ligand ' + sam + '_het.pdb -extra_res_fa pep.params -lig_grid -pocket_static_grid true -round_pocketGrid_center false -multiple_origin' 
#	print cmd4
	subprocess.call(cmd4,shell=True)
	return sips

def prep_lig(inp):
	cmd = GEAR + '/reduce -Trim ' + inp + '.pdb > ' + inp + '_red.pdb'
	subprocess.call(cmd,shell=True)	
	cmd2 = 'babel -ipdb ' + inp + '_red.pdb -omol2 ' + inp + '_red.mol2'
	subprocess.call(cmd2,shell=True)
	cmd4 = GEAR + '/balloon -f ' + LIB + '/MMFF94.mff --nconfs 100 --nGenerations 300 --rebuildGeometry ' + inp + '_pep_red.mol2 conformer.mol2 --noGA'
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

def reduce(pdbx,xch):
	cmd = '''python %s/aareduce.py %s %s-aa.pdb --chain %s --dumppatch --pdb2pqr > receptor.mapping
python %s/aareduce.py %s-aa.pdb %s-heavy.pdb --heavy --chain %s --readpatch  > /dev/null
python %s/reduce.py %s-aa.pdb %sr.pdb --chain %s > /dev/null'''%(ATTRACT_ALLATOM,pdbx,pdbx.split('.')[0],xch,ATTRACT_ALLATOM,pdbx.split('.')[0],pdbx.split('.')[0],xch,ATTRACT_TOOL,pdbx.split('.')[0],pdbx.split('.')[0],xch)
#	print cmd
	subprocess.call(cmd,shell=True)

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
	sam1 = sys.argv[2] # pep_seq
	sam2 = sys.argv[3] # N.of.CPU

	ROSETTA_BIN = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/source/bin'
	ROSETTA_DB = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/database'
	ATTRACT = '/awork06-1/attract'
	ATTRACT_BIN = '/awork06-1/attract/bin'
	ATTRACT_TOOL = '/awork06-1/attract/tools'
	ATTRACT_ALLATOM = '/awork06-1/attract/allatom'
	GEAR = '/awork06-1/neoscan_gear'

	wdir = os.getcwd()

	try:
		if not os.path.exists(sam + '_' + sam1):
			os.mkdir(sam + '_' + sam1,0777)
	except OSError:
		pass

	try:
		if not os.path.exists(sam + '_' + sam1 + '/result'):
			os.mkdir(sam + '_' + sam1 + '/result')
	except OSError:
		pass

	RES = 'result'

	# parameter
	params = ATTRACT + '/attract.par ' + sam + '_recr.pdb model/peptide1r.pdb --fix-receptor --ens 2 ens.list --gravity 2'
	scoreparams = ATTRACT + '/attract.par ' + sam + '_recr.pdb model/peptide1r.pdb --fix-receptor --ens 2 ens.list'
	gridparams = ' --grid 1 receptorgrid.gridheader'
	parals = '--np ' + sam2 + ' --chunks ' + sam2 

	shutil.copy(sam + '.pdb',sam + '_' + sam1)
	os.chdir(sam + '_' + sam1)
	ch = []
	rch = ''
	pch = ''
	with open(sam + '.pdb','r') as ffx:
		lines = ffx.readlines()
		for line in lines:
			if line.startswith('ATOM') > 0 and line[12:16].find('CA') > 0:
				ch.append(line[21])
	results = Counter(ch)
	max_ch = 0
	k = 0
	for key in results:
		if k == 0 :
			max_ch = results[key]
			rch  = key
		else:
			if results[key] > max_ch :
				rch = key
				pch = key1
			else:
				pch = key
		key1 = key
		k = k + 1
		
	prep_rec(sam,rch,pch)
	pep_build(sam1)
	reduce(sam + '_rec.pdb',rch)
	reduce('model/peptide1.pdb',pch)
	reduce('model/peptide2.pdb',pch)
	reduce('model/peptide3.pdb',pch)

	# For reference PDB( receptor, peptide)
	cmd = 'python ' + ATTRACT_ALLATOM + '/aareduce.py ' + sam + '_rec.pdb ' + sam + '_rec_ref.pdb --heavy --pdb2pqr > /dev/null'
	subprocess.call(cmd,shell=True)
	cmd1 = 'python ' + ATTRACT_ALLATOM + '/aareduce.py ' + sam + '_pep.pdb ' + sam + '_pep_ref.pdb --heavy --pdb2pqr > /dev/null'
	subprocess.call(cmd1,shell=True)

	# Generate starting structures & ensemble search
	cmd2 = 'python ' + ATTRACT_TOOL + '/randsearch.py 2 100000 --fix-receptor > randsearch.dat'
	subprocess.call(cmd2,shell=True)
	cmd3 = 'python ' + ATTRACT_TOOL + '/ensemblize.py randsearch.dat 3 2 all  > randsearch-ens2.dat'
	subprocess.call(cmd3,shell=True)

	grids = []
	with open('model/peptide1r.pdb','r') as ff:
		lines = ff.readlines()
		for line in lines:
			grids.append(line[57:59])

	grids = remove_duplicates(grids)
	grids.sort()
	print grids
	with open('receptorgrid.alphabet','w') as rf:
		for grid in grids:
			text = '%s\n'%(grid)
			rf.write(text)

	cmd4 = ATTRACT_BIN + '/make-grid-omp ' + sam + '_recr.pdb ' + ATTRACT + '/attract.par 5.0 7.0 receptorgrid.gridheader --shm --alphabet receptorgrid.alphabet'
	print cmd4
	subprocess.call(cmd4,shell=True)

	cmd5 = 'python ' + ATTRACT + '/protocols/attract.py randsearch-ens2.dat ' + params + ' ' + gridparams + ' --vmax 1000 ' + parals + ' --output out_' + sam + '_' + sam1 + '.dat'
	subprocess.call(cmd5,shell=True)

	cmd6 = 'python ' + ATTRACT + '/protocols/attract.py out_' + sam + '_' + sam1 + '.dat ' + scoreparams + ' --rcut 50.0 ' + parals + ' --output out_' + sam + '_' + sam1 + '.score'
	subprocess.call(cmd6,shell=True)

	cmd7 = 'python ' + ATTRACT_TOOL + '/fill-energies.py out_' + sam + '_' + sam1 + '.dat out_' + sam + '_' + sam1 + '.score > out_' + sam + '_' + sam1 + '-scored.dat'
	subprocess.call(cmd7,shell=True)

	cmd8 = 'python ' + ATTRACT_TOOL + '/sort.py out_' + sam + '_' + sam1 + '-scored.dat > out_' + sam + '_' + sam1 + '-sorted.dat'
	subprocess.call(cmd8,shell=True)

	cmd9 = ATTRACT_BIN + '/fix_receptor out_' + sam + '_' + sam1 + '-sorted.dat 2 --ens 0 3 | python ' + ATTRACT_TOOL + '/fill.py /dev/stdin out_' + sam + '_' + sam1 + '-sorted.dat > out_' + sam + '_' + sam1 + '-sorted.dat-fixre'
	subprocess.call(cmd9,shell=True)

	cmd10 = ATTRACT_BIN + '/deredundant out_' + sam + '_' + sam1 + '-sorted.dat-fixre 2 --ens 0 3 | python ' + ATTRACT_TOOL + '/fill-deredundant.py /dev/stdin out_' + sam + '_' + sam1 + '-sorted.dat-fixre > out_' + sam + '_' + sam1 + '-sorted-dr.dat'
	subprocess.call(cmd10,shell=True)

	cmd11 = ATTRACT_TOOL + '/top out_' + sam + '_' + sam1 + '-sorted-dr.dat 50 > out_' + sam + '_' + sam1 + '-pre-iattract.dat'
	subprocess.call(cmd11,shell=True)

	cmd12 = 'python ' + ATTRACT + '/protocols/iattract.py --infinite out_' + sam + '_' + sam1 + '-pre-iattract.dat ' + ATTRACT_ALLATOM +  '/allatom.par ' + sam +  '_rec-aa.pdb model/peptide1-aa.pdb --cdie --epsilon 10 --fix-receptor --icut 5.0 --np 16 --ens 2 ens-aa.list --name iattract-' + sam + '_' + sam1 + ' --output out_' + sam + '_' + sam1 + '-iattract.dat'
	subprocess.call(cmd12,shell=True)

	cmd13 = ATTRACT_TOOL + '/top out_' + sam + '_' + sam1 + '-iattract.dat 50 > out_' + sam + '_' + sam1 + '-top50.dat'
	subprocess.call(cmd13,shell=True)

	cmd14 = ATTRACT_BIN + '/collect out_' + sam + '_' + sam1 + '-top50.dat ' + sam + '_rec-aa.pdb model/peptide1-aa.pdb --name iattract-' + sam + '_' + sam1 + ' --ens 2 ens-aa.list  > out_' + sam + '_' + sam1 + '-top50.pdb'
	subprocess.call(cmd14,shell=True)
	
	idx = 1
	with open('out_' + sam + '_' + sam1 + '-top50.pdb','r') as ff:
		lines = ff.readlines()
		for line in lines:
			if line.startswith('MODEL') > 0 :
				pdb_out = open(RES + '/model_' + str(idx) + '.pdb','w')
				print >> pdb_out,'%s'%(line[:-1])
			elif line.startswith('ENDMDL') > 0 :
				print >> pdb_out,'END'
				pdb_out.close()
				idx = idx + 1
			else:
				print >> pdb_out,'%s'%(line[:-1])
	
	pdbs = glob.glob('result/*.pdb')
	for pdb in pdbs:
		cmd = GEAR + '/pep_rmsd ' + sam + '.pdb ' + pdb + ' bb ' + pch + ' >> rmsd.txt'
		subprocess.call(cmd,shell=True)

	with open('result_stat.txt','w') as rr:
		rr.write('PDB\tSEED\tEnergy\n')
		with open('out_' + sam + '_' + sam1 + '-top50.dat','r') as xx:
			lines = xx.readlines()
			for line in lines:
				if line.find('pivot') < 0 and line.find('centered') < 0 and line.find('Command') < 0 :
					if line.startswith('###') > 0 and line.find('SEED') > 0 :
						cols = line[:-1].split(' ') 
						rr.write('\t%s'%(cols[2]))
					elif line.startswith('##') > 0 and line.find('Energy') > 0 :
						cols = line[:-1].split(' ')
						rr.write('\t%s'%(cols[2]))
						rr.write('\n')
					elif line.startswith('#') > 0 :
						if line.startswith('##') > 0 :
							continue
						rr.write('model_%s'%(line[:-1][1:len(line[:-1])]))

	df = pd.read_csv('rmsd.txt',sep='\t',names=('REF','PDB','RMSD'))
	df1 = pd.read_csv('result_stat.txt',sep='\t')
	df2 = pd.merge(df,df1)
	df2.sort_values(['RMSD'],ascending=[True]).to_csv('final_result_stat.txt',sep='\t',index=False,na_rep='-')

	os.chdir(wdir)
					
#	cmd15 = 'python ' + ATTRACT_BIN + '/lrmsd.py out_' + sam + '_' + sam1 + '-iattract.dat model/peptide1-aa.pdb ' + sam + '_pep_ref.pdb --name iattract-' + sam + '_' + sam1 + ' --ens 2 ens-aa.list --receptor ' + sam + '_rec-aa.pdb > out_' + sam + '_' + sam1 + '-iattract.lrmsd' 
#	subprocess.call(cmd15,shell=True)
