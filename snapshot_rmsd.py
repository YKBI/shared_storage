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
import getInterfaceResidues
import Bio.PDB
from operator import itemgetter

def remove_duplicates(li):
	my_set = set()
	res = []
	for e in li:
		if e not in my_set:
			res.append(e)
			my_set.add(e)
	return res

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

pep_bind = []
if os.path.exists(feature):
	with open(feature,'r') as pf:
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

cwd = os.getcwd()
os.chdir(cwd)

nres = 0
resid = -99999
try:
	initial_pdb = open(glob.glob('%s/prep/*_initial_solv.pdb' % (cwd))[0], 'r')
except:
	org_initial_pdb = glob.glob('%s/prep/*_tleap.pdb' % (cwd))[0]
	shutil.copy(org_initial_pdb, '%s/prep/%s_initial_solv.pdb' % (cwd, tag))
	initial_pdb = open('cwd/prep/%s_initial_solv.pdb' % (cwd, tag), 'r')

atoms = initial_pdb.readlines()
for at in atoms:
	if at.startswith('ATOM'):
	#	print at
		this_resid = int(at[22:26])
		this_resnam = at[17:20]
		if this_resnam not in ['WAT', 'Na+', 'Cl-']:
			if this_resid != resid:
				nres += 1
				resid = this_resid
			else:pass
		else:
			break
str_comp_resid = 1
end_comp_resid = nres
str_prot_resid = 1
end_prot_resid = nres - plen
if plen > 0:
	lig_resid = []
	for i in range(1, plen+1):
		lig_resid.append(nres - i + 1)
	lig_resid.sort()
rec_str_resid = 1
rec_end_resid = end_prot_resid
if plen > 0 :
	lig_str_resid = lig_resid[0]
	lig_end_resid = lig_resid[-1]

traj_dirs = glob.glob('traj_*')
traj_dirs.sort()
if stTrajId+len(traj_dirs)-1 < numTraj:
	edTrajId = len(traj_dirs)
else:
	edTrajId = stTrajId + numTraj - 1

cpu_num = multiprocessing.cpu_count()

for i in range(stTrajId-1, edTrajId):
	this_traj = traj_dirs[i]
	os.chdir('%s/%s' % (cwd,this_traj))
	production = 'production'
	cwdt = os.getcwd()
#	print cwdt
	## trajectory list
	os.chdir('%s/production'%(cwdt))
	crd_files = glob.glob('%s_*.crd' % (tag))
	crd_files.sort()
	crd_loc = production
	## generate each snapshot
	os.chdir('%s/production'%(cwdt))
	first_ref = glob.glob('%s_*_md0000.rst' % (tag))[0]
	with open('md.in','r') as mdinfile:
		mdopts = mdinfile.readlines()
		for opt in mdopts:
			if opt.split()[0] == 'dt':
				dt = float(opt.split('=')[-1].split(',')[0])
			elif opt.split()[0] == 'nstlim':
				nstlim = int(opt.split('=')[-1].split(',')[0])
			elif opt.split()[0] == 'ntwx':
				ntwx = int(opt.split('=')[-1].split(',')[0])
			else: pass
	crd_t = float(dt*nstlim)
	crd_1 = float(dt*ntwx)
	freq = int(crd_t/crd_1)
	
	os.chdir('%s' % (cwdt))
	if (not os.path.exists(pdb_from_prod)):
		os.mkdir(pdb_from_prod)
		os.chdir(pdb_from_prod)

		anal_in = open('gen_pdb.in','w')
		print >> anal_in, '# set topology'
		print >> anal_in, 'parm ../%s_solv.top [%s_top]' % ( tag, tag)
		print >> anal_in, '# reference structure 1'
		print >> anal_in, 'reference ../../%s_initial_solv.crd parm [%s_top] [XRAY]' % (tag,tag)
		print >> anal_in, '# reference structure 2'
		print >> anal_in, 'reference ../production/%s parm [%s_top] [FIRST]' % ( first_ref, tag)

		idx = 0 
		for cf in crd_files:
			if idx == numcrd :
				break
			print >> anal_in, 'trajin ../%s/%s 1 %d parm [%s_top]' % (crd_loc,cf,freq,tag)
			idx = idx + 1

		print >> anal_in, 'autoimage'
		print >> anal_in, 'strip :WAT'
		print >> anal_in, 'strip :Na+'
		print >> anal_in, 'strip :Cl-'
		print >> anal_in, 'trajout %s.%s.pdb pdb multi' % (tag, this_traj)
		print >> anal_in, 'rms BB_rms @N,CA,C reference out rmsd.dat'
		print >> anal_in, 'rms Rec_rms :%s-%s@N,CA,C reference out rec_rmsd.dat' % ( rec_str_resid,rec_end_resid)
		print >> anal_in, 'atomicfluct out full_rmsf.dat byres bfactor'
		if len(pep_bind) > 0 :
			this_pep_bind_string1 = ','.join([str(x) for x in pep_bind])
			print >> anal_in, 'rms pep_bind_rms :%s@N,CA,C reference out pep_bind_rmsd.dat' % ( this_pep_bind_string1)
		if plen > 0 :
			if len(pep_bind) > 0 :
				print >> anal_in, 'rms pep_bind_rms_lig :%s,%s-%s reference out pep_bind_rmsd_lig.dat' % ( this_pep_bind_string1,lig_str_resid,lig_end_resid)
			print >> anal_in, 'rms Lig_rms :%s-%s reference out lig_rmsd.dat' % ( lig_str_resid,lig_end_resid)
			print >> anal_in, 'atomicfluct :%s-%s out lig_rmsf.dat bymask bfactor' %( lig_str_resid,lig_end_resid)
		anal_in.close()

		run_sh = open('run_gen_pdb.sh','w')
		print >> run_sh, ('''#!/bin/sh
		cpptraj -i gen_pdb.in''')
		run_sh.close()
		st = os.stat('run_gen_pdb.sh')
		os.chmod('run_gen_pdb.sh', st.st_mode | stat.S_IEXEC)
		subprocess.call(['./run_gen_pdb.sh'])

#		trms = []
#		with open('rmsd.dat','r') as rms :
#			for line in rms :
#				if (not line.startswith('#')):
#					rmsx = line.split()
#					trms.append(rmsx[0:])
#
#		trms1 = []
#		for trmsx in trms :
#			trmsx = map(float,trmsx)
#			trms1.append(trmsx)
#
#		trms1.sort(key=lambda x:x[1])
	else :
		os.chdir(pdb_from_prod)
#		trms = []
#		with open('rmsd.dat','r') as rms :
#			for line in rms :
#				if (not line.startswith('#')):
#					rmsx = line.split()
#					trms.append(rmsx[0:])
#
#		trms1 = []
#		for trmsx in trms :
#			trmsx = map(float,trmsx)
#			trms1.append(trmsx)
#
#		trms1.sort(key=lambda x:x[1])

	os.chdir('%s' % (cwdt))

	pdb_list = glob.glob('%s/*.pdb.*' % (pdb_from_prod))
	pdb_list.sort()

	rpdb_ids = []
	for pdb in pdb_list :
		rpdb_ids.append(int(pdb.split('.')[3]))

	rpdb_ids.sort()
#	print rpdb_ids
#	print len(rpdb_ids)

#	rpdb_ids =[]
#	kk = 0
#	for trmsx1 in trms1 :
#		if kk == numsnap :
#			break
#		rpdb_ids.append(trmsx1[0])	
#		kk += 1

#	rpdb_ids = map(int, rpdb_ids)
#	rpdb_ids.sort()

#	print rpdb_ids
#	print len(rpdb_ids)

	if (not os.path.exists(ML_feature_dir)):
		os.mkdir(ML_feature_dir)
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

	print "Extracting interacting residues from MD trajectory - %s" % this_traj
	conf_id=1
	for rpdb_id in rpdb_ids :
		rpdb = "../pdb_from_prod/" + tag + "." + this_traj + ".pdb." + str(rpdb_id)
		interacting_residue_pair = getInterfaceResidues.getInteractingResidueWithinRange(rpdb, rec_str_resid, rec_end_resid,lig_str_resid, lig_end_resid, dist_cut)
		tmp_rec = []
		tmp_lig = []
		tmp_chg_rec = []
		tmp_plr_rec = []
		tmp_npl_rec = []
		for res_pair in interacting_residue_pair:
			if (for_cnn):
				tmp_rec.append(res_pair[0][3:])
				if res_pair[0][:3] in polar_residue:
					tmp_plr_rec.append(res_pair[0][3:])
				elif res_pair[0][:3] in charged_residue:
					tmp_chg_rec.append(res_pair[0][3:])
				elif (res_pair[0][:3] not in charged_residue) and (res_pair[0][:3] not in polar_residue):
					tmp_npl_rec.append(res_pair[0][3:])
				else: pass
				if (res_pair[1][3:] not in tmp_lig):
					tmp_lig.append(res_pair[1][3:])
			if int(res_pair[0][3:]) not in rec_resid_all:
				rec_resid_all.append(int(res_pair[0][3:]))
			if int(res_pair[1][3:]) not in lig_resid_all:
				lig_resid_all.append(int(res_pair[1][3:]))

			if (res_pair[0][:3] in polar_residue) and (int(res_pair[0][3:]) not in polar_rec_resid_all):
				polar_rec_resid_all.append(int(res_pair[0][3:]))
			elif (res_pair[0][:3] in charged_residue) and (int(res_pair[0][3:]) not in charged_rec_resid_all):
				charged_rec_resid_all.append(int(res_pair[0][3:]))
			elif (res_pair[0][:3] not in charged_residue) and (res_pair[0][:3] not in polar_residue) and (int(res_pair[0][3:]) not in nonpolar_rec_resid_all):
				nonpolar_rec_resid_all.append(int(res_pair[0][3:]))
			else:pass

		if (for_cnn):
			lig_resid_for_conf[rpdb_id]            = remove_duplicates(tmp_lig)
			rec_resid_for_conf[rpdb_id]            = remove_duplicates(tmp_rec)
			charged_rec_resid_for_conf[rpdb_id]    = remove_duplicates(tmp_chg_rec)
			polar_rec_resid_for_conf[rpdb_id]      = remove_duplicates(tmp_plr_rec)
			nonpolar_rec_resid_for_conf[rpdb_id]   = remove_duplicates(tmp_npl_rec)
			lig_resid_for_conf[rpdb_id].sort()
			rec_resid_for_conf[rpdb_id].sort()
			charged_rec_resid_for_conf[rpdb_id].sort()
			polar_rec_resid_for_conf[rpdb_id].sort()
			nonpolar_rec_resid_for_conf[rpdb_id].sort()
		#	print "%d\t%s\t%s\t%s\t%s\t%s"%(rpdb_id,lig_resid_for_conf[rpdb_id],rec_resid_for_conf[rpdb_id],charged_rec_resid_for_conf[rpdb_id],polar_rec_resid_for_conf[rpdb_id],nonpolar_rec_resid_for_conf[rpdb_id])

	rec_resid_all.sort()
	lig_resid_all.sort()
	polar_rec_resid_all.sort()
	charged_rec_resid_all.sort()
	nonpolar_rec_resid_all.sort()
	rec_resid_string = ','.join([str(x) for x in rec_resid_all])
	lig_resid_string = ','.join([str(x) for x in lig_resid_all])
	polar_rec_resid_string = ','.join([str(x) for x in polar_rec_resid_all])
	charged_rec_resid_string = ','.join([str(x) for x in charged_rec_resid_all])
	nonpolar_rec_resid_string = ','.join([str(x) for x in nonpolar_rec_resid_all])

	run_feature_sh = open('run_feature.sh', 'w')
	print >> run_feature_sh, '#!/bin/sh\n'
	pdb_id = 1
	idx=0;
	print freq
	for i in range(len(crd_files)):
		if i == numcrd :
			break
		for j in range(1, freq+1):
	#		print "X11\t%d\t%d\t%d"%(idx,j,pdb_id)
	#		print "X1\t%d\t%d\t%d\t%d"%(idx,j,rpdb_ids[idx],pdb_id)
			if rpdb_ids[idx] == pdb_id :
		#		print "X2\t%d\t%d\t%d\t%d"%(idx,j,rpdb_ids[idx],pdb_id)
				feature_b_in = open('feature_%05d_%04d_bfactor.in' % (i+1, j), 'w')
				print >> feature_b_in, '# set topology'
				print >> feature_b_in, 'parm ../%s_solv.top [%s_top]' % ( tag, tag)
				print >> feature_b_in, '# reference structure 1'
				print >> feature_b_in, 'reference ../../%s_initial_solv.crd parm [%s_top] [XRAY]' % (tag,tag)
				print >> feature_b_in, '# reference structure 2'	
				print >> feature_b_in, 'reference ../production/%s parm [%s_top] [FIRST]' % ( first_ref, tag)
				print >> feature_b_in, 'trajin ../%s/%s 1 %d 1 parm [%s_top]' % (crd_loc, crd_files[i], j, tag)
				print >> feature_b_in, 'strip :WAT'
				print >> feature_b_in, 'strip :Na+'
				print >> feature_b_in, 'strip :Cl-'
				this_rec_resid_string = ','.join([str(x) for x in rec_resid_for_conf[pdb_id]])				
				this_lig_resid_string = ','.join([str(x) for x in lig_resid_for_conf[pdb_id]])
				print >> feature_b_in, 'atomicfluct :%s out rec_rmsf_bymask_bfactor_%05d_%04d.dat bymask bfactor' % (this_rec_resid_string, i+1, j)
				print >> feature_b_in, 'atomicfluct :%s out lig_rmsf_bymask_bfactor_%05d_%04d.dat bymask bfactor' % (this_lig_resid_string, i+1, j)
				feature_b_in.close()

				feature_in = open('feature_%05d_%04d.in' % (i+1, j), 'w')
				print >> feature_in, '# set topology'
				print >> feature_in, 'parm ../%s_solv.top [%s_top]' % ( tag, tag)
				print >> feature_in, '# reference structure 1'
				print >> feature_in, 'reference ../../%s_initial_solv.crd parm [%s_top] [XRAY]' % (tag,tag)
				print >> feature_in, '# reference structure 2'
				print >> feature_in, 'reference ../production/%s parm [%s_top] [FIRST]' % ( first_ref, tag)
				print >> feature_in, 'trajin ../%s/%s %d %d 1 parm [%s_top]' % (crd_loc, crd_files[i], j, j, tag)
				print >> feature_in, 'strip :WAT'
				print >> feature_in, 'strip :Na+'
				print >> feature_in, 'strip :Cl-'
				print >> feature_in, 'energy :%d-%d out energy_all_%05d_%04d.dat bond angle dihedral nb14 nonbond' % (str_prot_resid, nres, i+1, j)
				print >> feature_in, 'energy :%s,%s out energy_rec_lig_%05d_%04d.dat bond angle dihedral nb14 nonbond' % (this_rec_resid_string, this_lig_resid_string, i+1, j)
				print >> feature_in, 'lie :%s :%s out lie_lig_rec_%05d_%04d.dat' % (this_lig_resid_string, this_rec_resid_string,  i+1, j)
				print >> feature_in, 'lie :%s :%d-%d out lie_all_%05d_%04d.dat' % (this_lig_resid_string, rec_str_resid, rec_end_resid, i+1, j)
				if (len(charged_rec_resid_for_conf[pdb_id]) != 0):	
					this_charged_rec_resid_string = ','.join(charged_rec_resid_for_conf[pdb_id])
					print >> feature_in, 'lie :%s :%s out lie_charge_%05d_%04d.dat' % (this_lig_resid_string, this_charged_rec_resid_string, i+1,j)
				if (len(polar_rec_resid_for_conf[pdb_id]) != 0):
					this_polar_rec_resid_string = ','.join(polar_rec_resid_for_conf[pdb_id])
					print >> feature_in, 'lie :%s :%s out lie_polar_%05d_%04d.dat' % (this_lig_resid_string, this_polar_rec_resid_string, i+1,j)
				if (len(nonpolar_rec_resid_for_conf[pdb_id]) != 0):
					this_nonpolar_rec_resid_string = ','.join(nonpolar_rec_resid_for_conf[pdb_id])
					print >> feature_in, 'lie :%s :%s out lie_nonpolar_%05d_%04d.dat' % (this_lig_resid_string, this_nonpolar_rec_resid_string, i+1,j)
				feature_in.close()
				print >> run_feature_sh, 'cpptraj -i feature_%05d_%04d.in\n' % (i+1, j)
				print >> run_feature_sh, 'cpptraj -i feature_%05d_%04d_bfactor.in\n' % (i+1, j)
				idx += 1
			#	print idx
			#	print "X3\t%d\t%d\t%d\t%d"%(idx,j,rpdb_ids[idx],pdb_id)
			pdb_id += 1
	run_feature_sh.close()
	st = os.stat('run_feature.sh')
	os.chmod('run_feature.sh', st.st_mode | stat.S_IEXEC)
 	subprocess.call(['./run_feature.sh'])

	feature_out = open('CNN_feature.csv','w')
	feature = ''
	tmp = []
	pdb_id=1
	idx=0
	for i in range(len(crd_files)):
		if i == numcrd :
			break
		for j in range(1, freq+1):
	#		print "X1\t%d\t%d\t%d"%(idx,rpdb_ids[idx],pdb_id)
			if rpdb_ids[idx] == pdb_id :
			#	print "X2\t%d\t%d\t%d"%(idx,rpdb_ids[idx],pdb_id)
				if (os.path.exists('rec_rmsf_bymask_bfactor_%05d_%04d.dat' % (i+1, j))):
					with open('rec_rmsf_bymask_bfactor_%05d_%04d.dat' % (i+1, j), 'r') as rec_rmsf:
						for line in rec_rmsf:
							if (not line.startswith('#')):
								rmsf = line.split()
								tmp.extend(rmsf[1:])
				else :
					tmp.extend(['0.0'])
				if (os.path.exists('lig_rmsf_bymask_bfactor_%05d_%04d.dat' % (i+1, j))):
					with open('lig_rmsf_bymask_bfactor_%05d_%04d.dat' % (i+1, j), 'r') as lig_rmsf:
						for line in lig_rmsf:
							if (not line.startswith('#')):
								rmsf = line.split()
								tmp.extend(rmsf[1:])
				else :
					tmp.extend(['0.0'])
				if (os.path.exists('energy_all_%05d_%04d.dat' % (i+1, j))):
					with open('energy_all_%05d_%04d.dat' % (i+1, j), 'r') as energy_all:
						for line in energy_all:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0'])
				if (os.path.exists('energy_rec_lig_%05d_%04d.dat' % (i+1, j))):
					with open('energy_rec_lig_%05d_%04d.dat' % (i+1, j), 'r') as energy_rl:
						for line in energy_rl:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0'])
				if (os.path.exists('lie_all_%05d_%04d.dat' % (i+1, j))):
					with open('lie_all_%05d_%04d.dat' % (i+1, j), 'r') as lie_all:
						for line in lie_all:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0'])
				if (os.path.exists('lie_lig_rec_%05d_%04d.dat' % (i+1, j))):
					with open('lie_lig_rec_%05d_%04d.dat' % (i+1, j), 'r') as lie_lr:
						for line in lie_lr:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0'])
				if (os.path.exists('lie_charge_%05d_%04d.dat' % (i+1, j))):
					with open('lie_charge_%05d_%04d.dat' % (i+1, j), 'r') as lie_charge:
						for line in lie_charge:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0'])
				if (os.path.exists('lie_polar_%05d_%04d.dat' % (i+1,j))):
					with open('lie_polar_%05d_%04d.dat' % (i+1, j), 'r') as lie_polar:
						for line in lie_polar:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0'])
				if (os.path.exists('lie_nonpolar_%05d_%04d.dat' % (i+1,j))):
					with open('lie_nonpolar_%05d_%04d.dat' % (i+1, j), 'r') as lie_nonpolar:
						for line in lie_nonpolar:
							if (not line.startswith('#')):
								ene = line.split()
								tmp.extend(ene[1:])
				else :
					tmp.extend(['0.0','0.0'])
				idx += 1
			pdb_id += 1
	feature = ','.join(tmp)
	print >> feature_out, feature
	feature_out.close()
