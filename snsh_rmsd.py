#!/usr/bin/python

import os
import sys
import csv
import glob
import stat
import shutil
import numpy as np
import argparse
import subprocess
import multiprocessing
#import getInterfaceResidues
import Bio.PDB
from operator import itemgetter
cwd = os.getcwd()
def remove_duplicates(li):
	my_set = set()
	res = []
	for e in li:
		if e not in my_set:
			res.append(e)
			my_set.add(e)
	return res

def initial_pdb(x,y):# x is cwd, y is tag
	nres = 0
	resid = -99999
	try: I_pdb = open(glob.glob(x + '/prep/*_initial_solv.pdb')[0],'r')
	except:
		O_pdb = glob.glob(x + '/prep/*_tleap.pdb')[0]
		shutil.copy(O_pdb,x + '/prep/' + y + '_initial_solv.pdb')
		I_pdb = open(cwd + '/prep/' + y + '_initial_solv.pdb','r')

	for at in I_pdb.readlines():
		if at.startswith('ATOM'):
			print at
			this_resid,this_resnam = int(at[22:26]),at[17:20]
			if this_resnam not in ['WAT','Na+','Cl-']:
				if this_resid != resid:
					nres += 1
					resid = this_resid
				else: pass
			else: break
	I_pdb.close()
	return nres

def ligand_resid(x,y): #x is nres, y is plen
	str_prot_resid = 1
	end_prot_resid = x - y
	lig_resid = []
	if y > 0:
		for i in range(1,y+1):
			lig_resid.append(x - i + 1)
		lig_resid.sort()
	rec_str_resid = 1
	rec_end_resid = end_prot_resid
	if y > 0:
		lig_str_resid = lig_resid[0]
		lig_end_resid = lig_resid[-1]
	return str_prot_resid,lig_str_resid,lig_end_resid,rec_str_resid,rec_end_resid

def traj(x):
	this_traj = traj_dirs[x]
	os.chdir(cwd + '/' + this_traj + '/production')
	crd_files = sorted(glob.glob(tag + '_*.crd'))
	first_ref = glob.glob(tag + '_*_md0000.rst')[0]
	with open('md.in','r') as mdinfile:
		for opt in mdinfile.readlines():
			if opt.split()[0] == 'dt': dt = float(opt.split('=')[-1].split(',')[0])
			elif opt.split()[0] == 'nstlim': nstlim = int(opt.split('=')[-1].split(',')[0])
			elif opt.split()[0] == 'ntwx': ntwx = int(opt.split('=')[-1].split(',')[0])
	crd_t = int(dt*nstlim)
	crd_1 = int(dt*ntwx)
	freq = int(crd_t/crd_1)

	return crd_files,first_ref,freq,this_traj
'''x = tag, y = first_ref, z = this_traj, a = rec_str_resid, b = rec_end_resid,
c = lig_str_resid, d = lig_end_resid, files = crd_files,h = freq'''
def mk_gen_pdb(x,y,z,a,b,c,d,files,h):#x = tag, y = first_ref
	idx = 0
	gen_pdb = ['# set topology',\
			   'parm ../' + x +'_solv.top [' + x + '_top]',\
			   '# reference structre 1',\
			   'reference ../../' + x + '_initial_solv.crd parm [' + x + '_top] [XRAY]',\
			   '# reference structure 2',\
			   'reference ../production/' + y + ' parm [' + x + '_top] [FIRST]']
	for cf in files:
		if idx == numcrd:break
		gen_pdb.append('trajin ../production/' + cf + ' 1 ' + str(h) + ' parm [' + x + '_top]')
	this_pep_bind_string1 = ','.join([str(i) for i in pep_bind])
	gen_pdb_v1 = ['autoimage',\
				  'strip :WAT',\
				  'strip :Na+',\
				  'strip :Cl-',\
				  'trajout ' + x +'.' + z + '.pdb pdb multi',\
				  'rms BB_rms @N,CA,C reference out rmsd.dat',\
				  'rms Rec_rms :' + str(a) + '-' + str(b) + '@N,CAmC reference out rec_rmsd.dat',\
				  'atomicfluct out full_rmsf.dat byres bfactor',\
				  'rms pep_bind_rms :' + this_pep_bind_string1 + '@N,CA,C reference out pep_bind_rmsd.dat']
	if plen > 0:
		gen_pdb_v2 = ['rms pep_bind_rms_lig :' + this_pep_bind_string1 + ',' + str(c) + '-' + str(d) + ' reference out pep_bind_rmsd_lig.dat',\
					  'rms Lig_rms :' + str(c) + '-' + str(d) + ' reference out lig_rmsd.dat',\
					  'atomicfluct :' + str(c) + '-' + str(d) + ' out lig_rmsf.dat bymask bfactor']

	in_script = gen_pdb + gen_pdb_v1 + gen_pdb_v2

	with open('gen_pdb.in','w') as gin:
		for line in '\n'.join(in_script):
			gin.write(line)
	with open('run_gen_pdb.sh','w') as rsh:
		rsh.write('''#!/bin/sh\ncpptraj -i gen_pdb.in''')

#a=ref, b=tag, c=this_traj, d=feature, e=spdb1, f=spdb2
def sheba(a,b,c,d,e,f,g):
	cmd_1 = GEAR + '/sheba_01 -x ' + a + ' ' + e
	cmd_2 = GEAR + '/sheba_01 -t ' + b + '.' + c + '.' + str(g) + '.trf ' + e
	cmd_3 = GEAR + '/rmsd ' + d + ' ' + a + ' ' + f + " bb >> ext_bb_rmsd.dat"
	cmd_4 = GEAR + '/rmsd ' + d + ' ' + a + ' ' + f + " all >> ext_all_rmsd.dat"
	subprocess.call(cmd_1,shell=True)
	subprocess.call(cmd_2,shell=True)
	shutil.move(b + '.' + c + '.' + str(g) + '.pdb.pdb',f)
	subprocess.call(cmd_3,shell=True)
	subprocess.call(cmd_4,shell=True)

#x = tag,y=numcrd
def ref_pdb_feat(x,y):
	wdir = os.getcwd()
	print wdir
	os.chdir('traj_1/production')
	cmd = 'ambpdb -p ../' + x + '_solv.top -c ' + x + '_300K_1_md00' + str(y) + '.rst > ' + wdir +'/'+ x + '-ref.pdb'
	subprocess.call(cmd,shell=True)
	os.chdir(wdir)
	ref_line = []
	with open(wdir +'/'+ x + '-ref.pdb','r') as F:
		for line in F.readlines():
			if 'Na+' in line.strip():pass
			elif 'WAT' in line.strip():pass
			else: ref_line.append(line.strip())
	with open(wdir+'/' + x + '-ref.pdb.x','w') as N:
		for line in ref_line:
			N.write(line + '\n')
	shutil.move(wdir+'/' + x + '-ref.pdb.x',wdir+'/' + x + '-ref.pdb')
	feature_list = []
	with open(wdir+'/' + x+ '-ref.pdb','r') as rf:
		for line in rf.readlines():
			if line.startswith('ATOM' or 'TER'):feature_list.append(int(line[23:27].strip()))
	for i in sorted(list(set(feature_list))):
		print i
	with open(wdir+'/' + x+ '.list','w') as ft:
		for line in sorted(list(set(feature_list))):
			ft.write(str(line) + '\n')
	return wdir

print ('''
###############################################################################
#                                                                             #
#       Extracting structural & energetic features for Machine-learning       #
#                                                                             #
#       ** Running analysis options about trajectory files                    #
#                                                                             #
###############################################################################
''')

parser = argparse.ArgumentParser(description = "Extractor of structural & energetic features for ML")
parser.add_argument('-t',dest='tag', help='name tag used in this simulation', default=None)
parser.add_argument('-n', '--num-trajectory', dest='num_traj', type=int, help='Number of considered trajectories, default: 1', default = 1)
parser.add_argument('-n0', '--starting-trajectory', dest='st_traj', type=int, help='Starting trajectory ID. default: 1', default=1)
parser.add_argument('-nc', '--num-crd', dest='num_crd', type=int, help='Number of considered crd, default: 10', default = 10)
parser.add_argument('-c','--crd_time', metavar='[sim. time of one crd file (ps)]', dest='crd_t', type=int, help='Simulation time of one crd file. default: None', default=None)
parser.add_argument('-f','--freq', metavar='[time of one divied crd (ps)]', dest='crd_1', type=int, help='Time of one divided trajectory from one crd file, default: writing frequency in crd file, default: None', default=None)
parser.add_argument('-len', '--pep-len', dest='len', type=int, default=9, help='peptide length, default: 9')
parser.add_argument('-cut', '--distance-cut', dest='dist_cut', type=float, default=5.0, help='Distance cutoff between receptor and ligand, default: 5.0')
parser.add_argument('-nbcut', '--nonbonding-cut', dest='nb_cut', type=float, default=12.0, help='Cutoff for calculating nonbonding interaction (vdw, elec), default: 12.0')
parser.add_argument('-mut-pos', '--mutation-points', dest='mut_pos', nargs='+', type=str, default=None, help='Mutation positions, default: None')
parser.add_argument('-nopdb', '--do-not-gen-pdb', dest='no_pdb', action='store_true', help='Do not generate pdbs from production crds')
parser.add_argument('-nofeature', '--do-not-gen-feature', dest='no_feature', action='store_true', help='Do not feature files from production crds')
parser.add_argument('-nostat', '--do-not-perform-statistics', dest='no_stat', action='store_true', help='Do not perform statistics with generated features')
parser.add_argument('-for-cnn', '--gen-data-for-convolutional-neural-network', dest='for_cnn', action='store_true', help='Generate data for convolutional neural network platform (do not average data points')
#parser.add_argument('-ns', '--num-snapshot', dest='num_snap', type=int, help='Number of snapshot, default: 100', default = 100)
parser.add_argument('-feat','--feature', dest='feature', help='peptide binding interface')
parser.add_argument('-r',dest='ref', help='External ref PDB for rmsd calc', default=None)

args = parser.parse_args()
tag            = args.tag
numTraj        = args.num_traj
stTrajId       = args.st_traj
numcrd         = args.num_crd
crd_t          = args.crd_t
crd_1          = args.crd_1
plen           = args.len
dist_cut       = args.dist_cut
nb_cut         = args.nb_cut
mut_pos        = args.mut_pos
no_pdb         = args.no_pdb
no_feature     = args.no_feature
no_stat        = args.no_stat
for_cnn        = args.for_cnn
#numsnap        = args.num_snap
feature        = args.feature
ref			   = args.ref
inp_dir = tag.split('-')[1].split('.')[0]
GEAR = '/awork08/93_hong/NGS_ARS/Gear'
pep_bind = []
wdir = ref_pdb_feat(tag,numcrd)


with open(wdir+'/' + tag + '.list','r') as pf:
	lines = pf.readlines()
	for line in lines:
		pep_bind.append(line[:-1])

try:
	bin_dir = os.environ['MD_BIN']
except:
	print ('ERROR: define environment variable "MD_BIN"')
	sys.exit()
local_script_dir  = '%s/md_local_script' % (bin_dir)

pdb_from_prod     = 'pdb_from_prod'
ML_feature_dir    = 'ML_features'

charged_residue   = ['ARG', 'LYS', 'GLU', 'ASP']
polar_residue     = ['SER', 'THR', 'ASN', 'GLN']
try:
	if not os.path.exists('/home/user1/' + inp_dir):
		cp_dir = 'cp -r ' + inp_dir + ' /home/user1/' + inp_dir
		subprocess.call(cp_dir)
except OSError: pass


os.chdir(cwd)
traj_dirs = sorted(glob.glob('traj_*'))
if (stTrajId + len(traj_dirs) - 1) < numTraj:
	edTrajId = len(traj_dirs)
else: edTrajId = stTrajId + numTraj - 1

nres = initial_pdb(cwd,tag)
str_prot_resid,lig_str_resid,lig_end_resid,rec_str_resid,rec_end_resid = ligand_resid(nres,plen)
for i in range(stTrajId-1,edTrajId):
	crd_files,first_ref,freq,this_traj = traj(i)
	os.chdir('../')
	try:
		if not os.path.exists(ML_feature_dir):
			os.mkdir(ML_feature_dir)
	except OSError:pass
	try:
		if not os.path.exists(pdb_from_prod):
			os.mkdir(pdb_from_prod,0777)
	except OSError:pass
	os.chdir(pdb_from_prod)
	mk_gen_pdb(tag,first_ref,this_traj,rec_str_resid,rec_end_resid,lig_str_resid,lig_end_resid,crd_files,freq)

	st = os.stat('run_gen_pdb.sh')
	os.chmod('run_gen_pdb.sh',st.st_mode|stat.S_IEXEC)
	#make snapshot pdbs
	subprocess.call(['./run_gen_pdb.sh'])

	rpdb_ids = sorted([int(pdb.split('.')[3]) for pdb in glob.glob('*.pdb.*')])
	for rpdb in rpdb_ids:
		spdb = tag + '.' + this_traj + '.pdb.' + str(rpdb)
		spdb1  = tag + '.' + this_traj + '.' + str(rpdb) + '.pdb'
		spdb2 = tag + '.' + this_traj + '.' + str(rpdb) + '_tr.pdb'
		if not os.path.exists(spdb2):
			shutil.copy(spdb,spdb2)
			sheba(ref,tag,this_traj,wdir+'/'+tag+'.list',spdb1,spdb2,rpdb)

'''
	os.chdir(ML_feature_dir)
	rec_resid_all = []
	lig_resid_all = []
	charged_rec_resid_all = []
	polar_rec_resid_all = []
	nonpolar_rec_resid_all = []

	rec_resid_for_conf = {}
	charged_rec_resid_for_conf = {}
	polar_rec_resid_for_conf = {}
	nonpolar_rec_resid_for_conf = {}
	lig_resid_for_conf = {}

	conf_id = 1
	print "Extracting interacting residues from MD trajectory - " + this_traj
	for rpdb_id in rpdb_ids:
		rpdb = '../pdb_from_prod/' + tag +

'''
