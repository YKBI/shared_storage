#!/usr/bin/env python

import os
import glob
import shutil
import argparse
import stat
import random
import subprocess

import sys

def ext_amd_parm(bb):
	os.chdir(wdir)
	os.system('cat %s %s > tmp.pdb'%(init_rec,init_lig))
	natom = 0
	nres = 0
	with open('tmp.pdb','r') as f:
		lines = f.readlines()
		for line in lines:
			if line.startswith('ATOM') or line.startswith('HETATM'):
				natom = natom + 1
				if line[12:16].find('CA') > 0 :
					nres = nres + 1
	if nres == 0 :
		nres = nres + 1
	EthreshP = 0
	alphaP = 0
	EthreshD = 0
	alphaD = 0
	idx = 0
	print '%d\t%d'%(natom,nres)
	with open(bb,'r') as f1:
		lines = f1.readlines()
		for line in lines:
			if line.find(' A V E R A G E S   O V E R') > 0 :
				idx = idx + 1
			elif line.find('-----------------------------') > 0 :
				idx = 0
			if idx > 0 and line.find('EPtot') > 0 :
				texts = ' '.join(line[1:].split()).split(' ')
				EthreshP = float(texts[8]) + 0.16 * natom
			elif idx > 0 and line.find('DIHED') > 0 :
				texts = ' '.join(line[1:].split()).split(' ')
				EthreshD = float(texts[8]) + 4*nres
	print '%d\t%d'%(natom,nres)
	alphaP = 0.16 * natom
	alphaD = 4 * nres * 0.25
	os.chdir('traj_1/production/')
	return EthreshP,EthreshD,alphaP,alphaD

def equil_corr(equil,tag):
	nstep = 0
	ntwx = 0
	frame = 0
	with open(equil + '.out') as f :
		lines = f.readlines()
		for line in lines:
			if line.find('NSTEP') > 0 :
				cols = ' '.join(line.split())
				col = cols.split(' ')
				nstep = int(col[2])

	with open(equil + '.in') as f1 :
		lines = f1.readlines()
		for line in lines:
			if line.find('ntwx') > 0 :
				ntwx = int(line[:-1].lstrip().replace(',','').split(' ')[2])

	frame = int(nstep/ntwx)

	anal_in = open('ext_pdb.in','w')
	print >> anal_in, 'parm ../%s_solv.top [%s_top]' % ( tag, tag)
	print >> anal_in, 'trajin %s_%s.crd parm [%s_top]' % ( tag, equil, tag)
	print >> anal_in, 'trajout %s_%s.rst7 onlyframes %d'% ( tag, equil,frame)
	anal_in.close()

	run_sh = open('ext_pdb.sh','w')
	print >> run_sh, '#!/bin/csh'
	print >> run_sh, 'cpptraj -i ext_pdb.in'
	run_sh.close()
	st = os.stat('ext_pdb.sh')
	os.chmod('ext_pdb.sh', st.st_mode | stat.S_IEXEC)
	subprocess.call(['./ext_pdb.sh'])

	with open(equil + '_t.in','w') as eif:
		with open(equil + '.in','r') as eif1:
			lines = eif1.readlines()
			for line in lines:
				if line.find('irest') > 0 :
					eif.write(' irest = 0,\n')
				elif line.find('ntx') > 0 :
					eif.write(' ntx = 1,\n')
				else:
					eif.write(line)

def check_equil(equil):
	idx = 0
	tdir = os.getcwd()
	if os.path.exists(equil) :
		with open(equil) as eof :
			lines = eof.readlines()
			for line in lines:
				if line.find('TIMINGS') > 0 :
					idx = idx + 1
	else:
		idx = 0
	return idx

print ('''
#####################################################################################
#                                                                                   #
#       PREPARATION _ Explicit water molecular dynamics simulation using AMBER      #
#                       for protein-protein, protein-peptide complex                #
#       @ input : protein_receptor.pdb                                              #
#                 protein_ligand.pdb or peptide_ligand.pdb                          #
#                                                                                   #
#       ** Prepare H-deleted structure                                              #
#       ** Check protonation state of Histidine                                     #
#                                                                                   #
#####################################################################################
''')
#
# arguments
#
parser = argparse.ArgumentParser(description="Prepare simulation of THIS target")
parser.add_argument('-ir','--init_rec_str',metavar='[pdb initially processed, receptor]', dest='init_rec',required=True, help="initially processed receptor pdb structure located in working directory")
parser.add_argument('-il','--init_lig_str',metavar='[pdb initially processed, ligand]', dest='init_lig',required=True, help="initially processed ligand pdb structure located in working directory")
parser.add_argument('-t', '--tag_name',metavar='[name tag]',dest='tag',required=True, help='name tag used in generating top and initial crd file')
parser.add_argument('-n', '--num_traj',metavar='[num of trajectory]',dest='num_traj', type=int, help='number of generated trajectories, default: 6', default=6)
parser.add_argument('-n0','--traj0'   ,metavar='[starting trajectory id]',dest='start_traj', type=int, help='ID num of first trajectory', default=1)

parser.add_argument('-mpi', '--amber-mpi', dest='mpi_version', action='store_true', help='using pmemd.MPI ?')
parser.add_argument('-np', '--num-cpu', dest='num_cpu', type=int, help='number of CPU cores in MPI version', default=12)

parser.add_argument('-fix-gpu-id', '--fix-gpu-device-id', dest='fixed_gpu_id', type=int, help='Fix gpu device id (device id: 0,1,..), default: None', default=None)
parser.add_argument('-ngpu', '--num-gpu', dest='num_gpu', type=int, help='Number of GPU cards, default: None ==> 2', default=None)
parser.add_argument('-wb','--wat_box' ,metavar='[residual distance of water box, A]',dest='water_box_residual',type=float, help='distance from surface of protein to generate water box, default: 10.0 A', default=12.0)
parser.add_argument('-ec','--ewald-cut', metavar='[Ewald cutoff distance, A]', dest='ewald_cut', type=int, help='Ewald cut-off distance, default: 12 A', default=12)
parser.add_argument('-bfe','--set-for-binding-free-energy', dest='set_bfe', action='store_true', help='Flags for setting binding free energy')
parser.add_argument('-fix','--fix-pep_teriminal', dest='set_fix', action='store_true', help='Flags for fixing peptide terminal')

# restraint wt
parser.add_argument('-strong-restraint', dest='strong_restraint', type=float, default=5.0, help='strong restraint, default: 5.0')
parser.add_argument('-relaxed-restraint', dest='relaxed_restraint', type=float, default=2.0, help='relaxed restraint, default: 2.0')
parser.add_argument('-low-restraint', dest='low_restraint', type=float, default=1.0, help='low restraint, default: 1.0')
parser.add_argument('-minimal-restraint-1', dest='minimal_restraint_1', type=float, default=0.1, help='minimal restraint, default: 0.1')
parser.add_argument('-minimal-restraint-2', dest='minimal_restraint_2', type=float, default=0.5, help='minimal restraint, default: 0.5')

# minimization
parser.add_argument('-min1-maxcyc', dest='min1_maxcyc', type=int, default=3000, help='Num. of max cycle in 1st minimization, default: 3000 (steps)')
parser.add_argument('-min1-restraint-wt', dest='min1_rest_wt', type=float, default=100.0, help='Restraint force constant in 1st minimization (Restraints on protein & bound ligands), default: 100.0 (kcal/(mol*Angstrom^2))')
parser.add_argument('-min2-maxcyc', dest='min2_maxcyc', type=int, default=4000, help='Num. of max cycle in 2nd minimization, default: 4000 (steps)')

## simulation running options
parser.add_argument('-dt', metavar='[simulation time step, fs]', dest='dt', type=float, default=2, help='time step for simulation, defaulti: 2 fs')
parser.add_argument('-heat-time', metavar='[time for heating, ns]', dest='heat_time', type=float, default=1.0, help='simulation time for 1st equilibration (heating), default: 1 ns ')
parser.add_argument('-eq-restraint-wt', dest='eq_rest_wt', type=float, default=150.0, help='Restraint force constant in equilibration (heating), default: 150.0 (kcal/mol*Angstrom^2))')
parser.add_argument('-target-temp', metavar='[heating target temperature, K]', dest='tar_temp', type=float, default=300, help='target temperature for heating, default : 300 K')
parser.add_argument('-equil-time', metavar='[time for equilibration, ns]', dest='eq_time', type=float, default=10.0, help='simulation time for 2nd equilibration (equilibration), default: 10 ns ')
parser.add_argument('-pre-prod-time', metavar='[time for pre-MD, ns]', dest='eq_time1', type=float, default=0.0, help='simulation time for 2nd equilibration (equilibration), default: 0 ns ')
parser.add_argument('-prod-total-time', metavar='[total time for production, ns]', dest='prod_total_time', type=float, default=100.0, help='total time of simulation for production, default: 100 ns')
parser.add_argument('-prod-freq-time', metavar='[time for sub production run, ns]', dest='prod_freq_time', type=float, default=0.1, help='freq time of simulation for production, default: 0.1 ns')

## simulation result options
# EQ : heating, equilibration
parser.add_argument('-wrt-eq-out-freq', metavar='[writing frequency in .out file, EQ process, ns]', dest='wrt_eq_out_freq', type=float, default=0.01, help='freq time of writing EQ simulation output in .out file, default: 0.01 ns')
parser.add_argument('-wrt-eq-crd-freq', metavar='[writing frequency in .crd file, EQ process, ns]', dest='wrt_eq_crd_freq', type=float, default=0.01, help='freq time of writing .crd during running EQ simulation, default: 0.01 ns')
parser.add_argument('-wrt-eq-rst-freq', metavar='[writing freqeuncy in .rst file, EQ process, ratio]', dest='wrt_eq_rst_freq', type=float, default=1.0, help='frequecy of writing .rst during running EQ simulation, default: 1.0, at the end of simulation.')

# PROD
parser.add_argument('-wrt-pd-out-freq', metavar='[writing frequency in .out file, PD process, ns]', dest='wrt_pd_out_freq', type=float, default=0.01, help='freq time of writing PD simulation output in .out file, default: 0.01 ns')
parser.add_argument('-wrt-pd-crd-freq', metavar='[writing frequency in .crd file, PD process, ns]', dest='wrt_pd_crd_freq', type=float, default=0.01, help='freq time of writing .crd during running PD simulation, default: 0.01 ns')
parser.add_argument('-wrt-pd-rst-freq', metavar='[writing freqeuncy in .rst file, PD process, ratio]', dest='wrt_pd_rst_freq', type=float, default=1.0, help='frequecy of writing .rst during running PD simulation, default: 1.0, at the end of simulation.')
parser.add_argument('-amd',dest='iamd',type=float,default=3,help='Acceleration Molecular Dynamic simulation')

args = parser.parse_args()

init_rec    = args.init_rec
init_lig    = args.init_lig
tag         = args.tag
num_traj    = args.num_traj
traj0       = args.start_traj
mpi_version = args.mpi_version
num_cpu     = args.num_cpu

fixed_gpu_id   = args.fixed_gpu_id
num_gpu        = args.num_gpu
fix_gpu_device_id = False
max_gpu_dev_num = 4
if (fixed_gpu_id != None) and (max_gpu_dev_num > fixed_gpu_id > -1):
   if (num_gpu == None):
      num_gpu = 1
   elif num_gpu < 1:
      print (" ** ERROR: You assigned incorrect number of GPU device")
   else: # if num_gpu is assigned, use it
      pass
   fix_gpu_device_id = True
elif (fixed_gpu_id != None) and (fixed_gpu_id < 0):
   print (" ** ERROR: You assigned inappropriate GPU device id. Plz, Check it")
   sys.exit()
else:
   pass
if (num_gpu == None): num_gpu = 2      # default: num_gpu = 2

water_box_residual = args.water_box_residual
ewald_cut   = args.ewald_cut 
set_bfe     = args.set_bfe
set_fix     = args.set_fix

strong_restraint  = args.strong_restraint
relaxed_restraint = args.relaxed_restraint
low_restraint     = args.low_restraint
minimal_restraint_1  = args.minimal_restraint_1
minimal_restraint_2  = args.minimal_restraint_2

min1_maxcyc = args.min1_maxcyc
min1_rest_wt = args.min1_rest_wt    # kcal/(mol*Angstrom^2)
min2_maxcyc = args.min2_maxcyc

# Simulation running option
dt          = args.dt         # fs
heat_time   = args.heat_time  # ns
eq_rest_wt  = args.eq_rest_wt # kcal/(mol*Angstrom^2)
target_temp = args.tar_temp   # K
eq_time     = args.eq_time    # ns
eq_time1	=args.eq_time1 #ns
wrt_eq_out_freq      = args.wrt_eq_out_freq     # ns
wrt_eq_crd_freq      = args.wrt_eq_crd_freq     # ns
wrt_eq_rst_freq      = args.wrt_eq_rst_freq     # ratio, wrt_time/total_time

prod_total_time      = args.prod_total_time     # ns
prod_freq_time       = args.prod_freq_time      # ns
wrt_pd_out_freq      = args.wrt_pd_out_freq     # ns
wrt_pd_crd_freq      = args.wrt_pd_crd_freq     # ns
wrt_pd_rst_freq      = args.wrt_pd_rst_freq     # ratio, wrt_time/total_time
iamd                = args.iamd # acceleration md simulation
##########################################################
# inline function ??
def set_sim_options(filename, seed_number=None, process=None, str_res_id=None, end_res_id=None):
    shutil.move(filename, 'tmp.in')
    new_cont = []
    nstlim = None
    nptr = None
    ntwx = None
    ntwr = None
    fs2ns = 1000000     # femto sec to nano sec
    with open('tmp.in', 'r') as md_opt_file:
        md_opts = md_opt_file.readlines()
        for opt in md_opts:
            if ('SEED_NUMBER' in opt):
               new_cont.append(opt.replace('SEED_NUMBER','%d' % seed_number))
            elif ('EWALD_CUT' in opt):
               new_cont.append(opt.replace('EWALD_CUT', '%d' % ewald_cut))
            elif ('RESTRAINT_WT_MIN' in opt):
               new_cont.append(opt.replace('RESTRAINT_WT_MIN', '%.1f' % min1_rest_wt))
            elif ('MIN1_MAXCYC' in opt):
               new_cont.append(opt.replace('MIN1_MAXCYC', '%d' % min1_maxcyc))
            elif ('MIN1_NCYC' in opt):
               new_cont.append(opt.replace('MIN1_NCYC', '%d' % (min1_maxcyc/2)))
            elif ('MIN2_MAXCYC' in opt):
               new_cont.append(opt.replace('MIN2_MAXCYC', '%d' % min2_maxcyc))
            elif ('MIN2_NCYC' in opt):
               new_cont.append(opt.replace('MIN2_NCYC', '%d' % (min2_maxcyc/2)))
            elif ('TARGET_TEMP' in opt):
               new_cont.append(opt.replace('TARGET_TEMP', '%6.1f' % target_temp))
            elif ('TOTAL_HEAT_STEPS' in opt):
               nstlim = int(heat_time * fs2ns/dt)
               new_cont.append(opt.replace('TOTAL_HEAT_STEPS', '%d' % nstlim))
            elif ('TOTAL_EQ_STEPS' in opt):
               if process == 'PMD':
                  nstlim = int(eq_time1 * fs2ns/dt)
               else:
                  nstlim = int(eq_time * fs2ns/dt)
               new_cont.append(opt.replace('TOTAL_EQ_STEPS', '%d' % nstlim))
            elif ('TOTAL_PROD_STEPS' in opt):
               nstlim = int(prod_freq_time * fs2ns/dt)
               new_cont.append(opt.replace('TOTAL_PROD_STEPS', '%d' % nstlim))
            elif ('TIME_DT' in opt):
               _dt = 0.001*dt
               new_cont.append(opt.replace('TIME_DT', '%f' % _dt))
            elif ('STRONG_RESTRAINT' in opt):
               new_cont.append(opt.replace('STRONG_RESTRAINT', '%f' % strong_restraint))
            elif ('RELAXED_RESTRAINT' in opt):
               new_cont.append(opt.replace('RELAXED_RESTRAINT', '%f' % relaxed_restraint))
            elif ('LOW_RESTRAINT' in opt):
               new_cont.append(opt.replace('LOW_RESTRAINT', '%f' % low_restraint))
            elif ('MINIMAL_RESTRAINT_1' in opt):
               new_cont.append(opt.replace('MINIMAL_RESTRAINT_1', '%f' % minimal_restraint_1))
            elif ('MINIMAL_RESTRAINT_2' in opt):
               new_cont.append(opt.replace('MINIMAL_RESTRAINT_2', '%f' % minimal_restraint_2))
            elif ('RESTRAINT_WT_EQ' in opt):
               new_cont.append(opt.replace('RESTRAINT_WT_EQ', '%.1f' % eq_rest_wt))
            elif ('RESTRAINT_WT_PR' in opt):
               new_cont.append(opt.replace('RESTRAINT_WT_PR', '%.1f' % eq_rest_wt))
            elif ('WRT_OUT_FREQ' in opt):
               if process == 'HEAT' or process == 'EQ' or process == 'PMD':
                  nptr = int(wrt_eq_out_freq * fs2ns/dt)
               elif process == 'PROD':
                  nptr = int(wrt_pd_out_freq * fs2ns/dt)
               else:
                  print ("ERROR in setting out frequency")
                  sys.exit()
               new_cont.append(opt.replace('WRT_OUT_FREQ', '%d' % nptr))
            elif ('WRT_CRD_FREQ' in opt):
               if process == 'HEAT' or process == 'EQ' or process == 'PMD':
                  ntwx = int(wrt_eq_crd_freq * fs2ns/dt)
               elif process == 'PROD':
                  ntwx = int(wrt_pd_crd_freq * fs2ns/dt)
               else:
                  print ("ERROR in setting crd frequency")
                  sys.exit()
               new_cont.append(opt.replace('WRT_CRD_FREQ', '%d' % ntwx))
            elif ('WRT_RST_FREQ' in opt):
               if process == 'HEAT':
                  nstlim = int(heat_time * fs2ns/dt)
                  ntwr = int(wrt_eq_rst_freq*nstlim)
               elif process == 'EQ':
                  nstlim = int(eq_time * fs2ns/dt)
                  ntwr = int(wrt_eq_rst_freq*nstlim)
               elif process == 'PMD':
                  ntwr = int(wrt_eq_rst_freq*nstlim)
               elif process == 'PROD':
                  nstlim = int(prod_freq_time * fs2ns/dt)
                  ntwr = int(wrt_pd_rst_freq*nstlim)
		  if iamd != 0:
		      EthreshP, EthreshD, alphaP, alphaD = ext_amd_parm('traj_1/%s/equil_2.out'%equil_dir)
		      new_cont.append('  iamd=%d,\n'%iamd)
		      new_cont.append('  ethreshd=%f,alphad=%f,\n'%(EthreshD,alphaD))
		      new_cont.append('  ethreshp=%f,alphap=%f,\n'%(EthreshP,alphaP))
               else:
                  print ("ERROR in setting rst frequency")
                  sys.exit()
               new_cont.append(opt.replace('WRT_RST_FREQ', '%d' % ntwr))
            elif ('restraintmask' in opt):
               tmp_line = opt
               if ('STR_RES_NUM' in tmp_line):
                  tmp_line = tmp_line.replace('STR_RES_NUM', '%d' % str_res_id)
               if ('END_RES_NUM' in tmp_line):
                  tmp_line = tmp_line.replace('END_RES_NUM', '%d' % end_res_id)
               if ('HEAVY_ATOM_MASK' in tmp_line):
                  tmp_line = tmp_line.replace('HEAVY_ATOM_MASK', '*&!@H=')
               if ('BACKBONE_MASK' in tmp_line):
                  tmp_line = tmp_line.replace('BACKBONE_MASK', 'CA,C,O,N,H')
               new_cont.append(tmp_line)
            
            else: 
	       new_cont.append(opt)

      # re-write simulation parameter file
    with open(filename, 'w') as newf:
        for nc in new_cont:
            newf.write(nc)
    os.remove('tmp.in')

##########################################################

wdir = os.getcwd()  # working directory

prepdir = 'prep'
mindir  = 'min'
equil_dir = 'equilibration'
prod_dir  = 'production'
crd_dir   = 'crd'
pre_rsm_dir = 'pre_rsm'

#common = '/tools/bin_sim/common'
#bin_dir = '/tools/bin_sim'

try:
   common = '/tools/bin_sim/common' #os.environ['MD_COMMON']
except:
   print ('ERROR: define environment variable "MD_COMMON".')
   sys.exit()

try:
   bin_dir = '/tools/bin_sim'  #os.environ['MD_BIN']
except:
   print ('ERROR: define environment variable "MD_sim".')
   sys.exit()
water_type = 'TIP3PBOX'

gpu_dev_id = 0

if (set_bfe):     # if setting for binding free energy cal is required, set again
   num_traj = 30
   # Simulation running option
   heat_time = 0.1
   eq_time = 10      # long equilibration time maybe is required
   wrt_eq_rst_freq = 0.1   # ratio, wrt_rst_time / total_time

   prod_total_time = 10    # ns
   prod_freq_time  = 0.1   # ns
   wrt_pd_out_freq = 0.01  # ns
   wrt_pd_crd_freq = 0.01  # ns
   wrt_pd_rst_freq = 1     # ratio, wrt_rst_time / total_time
#
# prep dir
#
if (not os.path.exists(prepdir)):
    os.mkdir(prepdir)
else:
    pass
#---- prep dir
## prepparing ligand's frcmod and lib files and then 
## generating protein-ligand complex and its top, inpcrd files
os.chdir(prepdir)
shutil.copy('%s/%s'%(wdir,init_rec),'./%s'%(init_rec))
shutil.copy('%s/%s'%(wdir,init_lig),'./%s'%(init_lig))
ligname = init_lig.split('.')[0]       # ligand file name without file type
recname = init_rec.split('.')[0]       # receptor file name without file type
## merge receptor pdb and ligand pdb into complex pdb
os.system('cat %s.pdb >  %s.pdb' % (recname, tag))
os.system('cat %s.pdb >> %s.pdb' % (ligname, tag))
ff_type = 'leaprc.ff99SBildn'
amberhome = os.environ['AMBERHOME']
tleap_script = '''source %s/dat/leap/cmd/oldff/%s
complex = loadpdb %s.pdb
recept = loadpdb %s.pdb
ligand = loadpdb %s.pdb
''' % (amberhome, ff_type, tag, recname, ligname)

if (set_bfe):
   tleap_script += '''
set default PBRadii mbondi2
'''

tleap_script += '''
saveamberparm complex %s.prmtop %s.inpcrd
saveamberparm recept %s.prmtop %s.inpcrd
saveamberparm ligand %s.prmtop %s.inpcrd

solvateBox complex %s %6.1f iso
addIons complex Na+ 0
addIons complex Cl- 0
saveamberparm complex %s_solv.prmtop %s_solv.inpcrd
savepdb complex %s_initial_solv.pdb
quit
''' % (tag, tag, recname, recname, ligname, ligname,\
       water_type, water_box_residual,\
       tag, tag, tag)

tleap_script_f = open('tleap.script','w')
#print >> tleap_script_f, tleap_script
#print (tleap_script, end='', file=tleap_script_f)
tleap_script_f.write(tleap_script)
tleap_script_f.close()

tleap_run_script='''#!/bin/sh

tleap -f tleap.script -s
'''
tleap_run_script_f = open('run_tleap.sh','w')
#print >> tleap_run_script_f, tleap_run_script
#print (tleap_run_script, end='', file=tleap_run_script_f)
tleap_run_script_f.write(tleap_run_script)
tleap_run_script_f.close()
st = os.stat('run_tleap.sh')
os.chmod('run_tleap.sh', st.st_mode | stat.S_IEXEC)
subprocess.call('./run_tleap.sh',shell=True)

if (not os.path.lexists('%s_initial_solv.crd' % (tag))):
   os.symlink('%s_solv.inpcrd' % (tag), '%s_initial_solv.crd' % (tag))
if (not os.path.lexists('%s_solv.top')):
   os.symlink('%s_solv.prmtop' % (tag), '%s_solv.top' % (tag))

# get residue num of a receptor protein
pdbr = open('%s.pdb' % (recname),'r')
pdbrline = pdbr.readlines()
pdbr.close()
rstart = 1
rres = 0
resid = -99999
for p in pdbrline:
	if p.startswith('ATOM'):
		this_resid = int(p[22:26])
		if this_resid != resid:
			rres += 1
			resid = this_resid
		else:
			pass

# get residue num of a peptie
pdbl = open('%s.pdb' % (ligname),'r')
pdblline = pdbl.readlines()
pdbl.close()
lstart = 1
lres = 0
resid = -99999
for p in pdblline:
	if p.startswith('ATOM'):
		this_resid = int(p[22:26])
		if this_resid != resid:
			lres += 1
			resid = this_resid
		else:
			pass

lstart = lstart + rres
lres = lres + rres

# get residue num of a target protein, complex pdb
pdbf = open('%s.pdb' % (tag),'r')
pdbline = pdbf.readlines()
pdbf.close()
nres = 0
resid = -99999
for p in pdbline:
    if p.startswith('ATOM'):
        this_resid = int(p[22:26])
        if this_resid != resid:
            nres += 1
            resid = this_resid
        else:
            pass

#---- working directory
os.chdir(wdir)
# complex system in solvation
if (not os.path.lexists('./%s_initial_solv.crd' % (tag))):  # complex's solv crd 
    os.symlink('./%s/%s_initial_solv.crd'%(prepdir,tag),'./%s_initial_solv.crd'%(tag))
if (not os.path.lexists('./%s_solv.top' % (tag))):          # complex's solv top
    os.symlink('./%s/%s_solv.top'%(prepdir,tag), './%s_solv.top'%(tag))

# complex system in vaccum
if (not os.path.lexists('./%s.inpcrd' % (tag))):            # complex's solv crd 
    os.symlink('./%s/%s.inpcrd'%(prepdir,tag),'./%s.inpcrd' %(tag))
if (not os.path.lexists('./%s.prmtop' % (tag))):            # complex's solv top
    os.symlink('./%s/%s.prmtop'%(prepdir,tag), './%s.prmtop' %(tag))

# receptor system in vaccum
if (not os.path.lexists('./%s.inpcrd' % (recname))):        # receptor's vaccum crd
    os.symlink('./%s/%s.inpcrd'%(prepdir,recname), '%s.inpcrd'%(recname))
if (not os.path.lexists('./%s.prmtop' % (recname))):        # receptor's vaccum top
    os.symlink('./%s/%s.prmtop'%(prepdir,recname), '%s.prmtop'%(recname))

# ligand system in vaccum
init_lig_name = init_lig.split('.')[0]
if (not os.path.lexists('./%s.inpcrd' % (init_lig_name))):  # ligand's vaccum crd
    os.symlink('./%s/%s.inpcrd'%(prepdir,init_lig_name), '%s.inpcrd'%(init_lig_name))
if (not os.path.lexists('./%s.prmtop' % (init_lig_name))):  # ligand's vaccum top
    os.symlink('./%s/%s.prmtop'%(prepdir,init_lig_name), '%s.prmtop'%(init_lig_name))

#
# min dir
#
if (not os.path.exists(mindir)):
    os.mkdir(mindir)
else: pass
if (fix_gpu_device_id):
   gpu_dev_id = fixed_gpu_id
else:
   gpu_dev_id = 0
min_all_script = open('run_initialization.sh', 'w')
min_all_cont = '''#!/bin/sh
export CUDA_VISIBLE_DEVICES='%d'

cd min
./run_min1.sh
./run_min2.sh
cd ..
''' % (gpu_dev_id)
min_all_script.write(min_all_cont)
min_all_script.close()
st = os.stat('run_initialization.sh')
os.chmod('run_initialization.sh', st.st_mode | stat.S_IEXEC)

#---- min dir
os.chdir(mindir)
# prep for first minimization
shutil.copy('%s/simulation_in/min_1.in'%(common), './min_1.in')
set_sim_options('./min_1.in', process='MIN', str_res_id=1, end_res_id=nres)
#  # set residue num in first minimization input file
#  min1_f = open('temp_min_1.in','r')
#  min1 = min1_f.readlines()
#  min1_f.close()
#  os.remove('temp_min_1.in')
#  new_min1_f = open('min_1.in','w')
#  for i in range(len(min1)):
#      if min1[i].startswith('RES'):
#          min1[i] = 'RES 1 %d\n' % (nres)       # set residue range ( 1 ~ last residue )
#      else: pass
#      #print >> new_min1_f, min1[i][:-1]
#      #print (min1[i][:-1], end='\n', file=new_min1_f)
#      new_min1_f.write(min1[i][:-1])
#  new_min1_f.close()

if (mpi_version):
   min1_run_script = '''#!/bin/sh

mpirun -np %d pmemd.MPI -O -i min_1.in -o min_1.out -p ../%s_solv.top -c ../%s_initial_solv.crd -r %s_min1.rst -ref ../%s_initial_solv.crd 
''' % (num_cpu, tag, tag, tag, tag)
else:
   min1_run_script = '''#!/bin/sh

pmemd.cuda -O -i min_1.in -o min_1.out -p ../%s_solv.top -c ../%s_initial_solv.crd -r %s_min1.rst -ref ../%s_initial_solv.crd
''' % (tag, tag, tag, tag)
min1_run_f = open('run_min1.sh','w')
min1_run_f.write(min1_run_script)
min1_run_f.close()
st = os.stat('run_min1.sh')
os.chmod('run_min1.sh', st.st_mode | stat.S_IEXEC)

# prep for second minimization
shutil.copy('%s/simulation_in/min_2.in'%(common), './min_2.in')
set_sim_options('./min_2.in', process='MIN')

if (mpi_version):
   min2_run_script = '''#!/bin/sh

mpirun -np %d pmemd.MPI -O -i min_2.in -o min_2.out -p ../%s_solv.top -c %s_min1.rst -r %s_min2.rst
''' % (num_cpu, tag, tag, tag)
else:
   min2_run_script = '''#!/bin/sh

pmemd.cuda -O -i min_2.in -o min_2.out -p ../%s_solv.top -c %s_min1.rst -r %s_min2.rst
''' % (tag, tag, tag)
min2_run_f = open('run_min2.sh','w')
min2_run_f.write(min2_run_script)
min2_run_f.close()
st = os.stat('run_min2.sh')
os.chmod('run_min2.sh', st.st_mode | stat.S_IEXEC)
subprocess.call('./run_min1.sh',shell=True)
subprocess.call('./run_min2.sh',shell=True)

#---- working directory
os.chdir(wdir)
if (not os.path.exists('set_residue_info')):
    os.mkdir('set_residue_info')
if (not os.path.exists('set_atom_info')):
    os.mkdir('set_atom_info')
if (not os.path.exists('set_native_HP_core_contacts')):
    os.mkdir('set_native_HP_core_contacts')

#
#   each trajectory directory
#

###############################################
# gen seed number
seednums = []
while len(seednums) < num_traj:
    seed = random.randint(10000,100000)
    if seed not in seednums:
        seednums.append(seed)
seednums.sort()

###############################################
if (mpi_version):
   run_sim = open('run_sim_mpi.sh', 'w')
   run_sim.write('#!/bin/sh\n')
   for i in range(num_traj):
      if num_traj+1 < 10:
         traj_id = '%d' % (traj0+i)
      else:
         traj_id = '%02d' % (traj0+i)
      cont='''
cd ./traj_%s
cd ./equilibration
./run_equil.sh
cd ..
cd ./production
./run_production.sh
cd ../..

      ''' % (traj_id)
      run_sim.write(cont)
   run_sim.close()
   st = os.stat('run_sim_mpi.sh')
   os.chmod('run_sim_mpi.sh', st.st_mode | stat.S_IEXEC)
else: # gpu version
   # write run_sim.sh
   gpu_content = {}
   for i in range(num_traj):
       if num_traj+1 < 10:
          traj_id = '%d' % (traj0+i)
       else:
          traj_id = '%02d' % (traj0+i)
       cont = '''
cd ./traj_%s
cd ./equilibration
./run_equil.sh
cd ..
cd ./production
./run_production.sh
cd ../..

   ''' % (traj_id)

       if (not fix_gpu_device_id):
          if (i % num_gpu not in gpu_content):
             gpu_content[i % num_gpu] = [cont]
          else:
             gpu_content[i % num_gpu].append(cont)
       else:
          if (fixed_gpu_id not in gpu_content):
             gpu_content[fixed_gpu_id] = [cont]
          else:
             gpu_content[fixed_gpu_id].append(cont)

   for i in gpu_content.keys():
       if (not fix_gpu_device_id):
          run_sim_script = open('run_sim_gpu_%03d.sh' % i, 'w')
       else:
          if (not os.path.exists('run_sim_gpu_%03d.sh' % fixed_gpu_id)):
             run_sim_script = open('run_sim_gpu_%03d.sh' % fixed_gpu_id, 'w')
          else:
             os.remove('run_sim_gpu_%03d.sh' % (fixed_gpu_id))
             run_sim_script = open('run_sim_gpu_%03d.sh' % i, 'w')
       
       run_sim_script.write('#!/bin/sh\n\n')
       for c in gpu_content[i]:
          run_sim_script.write('%s\n' % (c))
       run_sim_script.close()
       if (not fix_gpu_device_id):
          st = os.stat('run_sim_gpu_%03d.sh' % i)
          os.chmod('run_sim_gpu_%03d.sh' % i, st.st_mode | stat.S_IEXEC)
       else:
          st = os.stat('run_sim_gpu_%03d.sh' % fixed_gpu_id)
          os.chmod('run_sim_gpu_%03d.sh' % fixed_gpu_id, st.st_mode | stat.S_IEXEC)

###############################################
# generate trajectory
for i in range(num_traj):
    
    ###############################################
    if num_traj+1 < 10:
       traj_id = '%d' % (traj0+i)
    else:
       traj_id = '%02d' % (traj0+i)
    this_traj = 'traj_%s' % (traj_id)
    if (not os.path.exists(this_traj)):
        os.mkdir(this_traj)
    
    os.chdir(this_traj)
    if (not os.path.exists(equil_dir))      : os.mkdir(equil_dir)
    if (not os.path.exists(prod_dir))       : os.mkdir(prod_dir)
    if (not os.path.exists(crd_dir))        : os.mkdir(crd_dir)
    if (not os.path.exists(pre_rsm_dir))    : os.mkdir(pre_rsm_dir)
    if (not os.path.lexists('./%s_solv.top'%(tag))):
        os.symlink('../%s_solv.top'%(tag), './%s_solv.top'%(tag))
    if (not os.path.lexists('./%s_min2.rst' % (tag))):
        os.symlink('../%s/%s_min2.rst' % (mindir,tag), './%s_min2.rst' % (tag))
    if (set_fix):
       os.symlink('../%s/%s_solv.inpcrd' % (prepdir,tag), './%s_initial_solv.crd' % (tag))
    
    ###############################################
    ##---- equilibration dir
    os.chdir(equil_dir)
    
    # modifying equil_1.in : heating
    shutil.copy('%s/simulation_in/tmp_equil_1.in' % (common), './equil_1.in')
    set_sim_options('./equil_1.in', seed_number=seednums[i], process='HEAT', str_res_id=1, end_res_id=nres)
    
    ### 2nd equilibration : equil.
    if (set_fix):    
       shutil.copy('%s/simulation_in/tmp_equil_2_fix.in' % (common), './equil_2.in')
       set_sim_options('./equil_2.in', process='EQ', str_res_id=lstart, end_res_id=lres)
       if iamd != 0 :
          shutil.copy('%s/simulation_in/tmp_equil_2_fix.in' % (common), './equil_3.in')
          set_sim_options('./equil_3.in', process='PMD', str_res_id=lstart, end_res_id=lres)
    else:
       shutil.copy('%s/simulation_in/tmp_equil_2.in' % (common), './equil_2.in')
       set_sim_options('./equil_2.in', process='EQ')
       if iamd != 0 :
          shutil.copy('%s/simulation_in/tmp_equil_2.in' % (common), './equil_3.in')
          set_sim_options('./equil_3.in', process='PMD')
       
    
    if (fix_gpu_device_id):
       gpu_dev_id = fixed_gpu_id
    else:
       gpu_dev_id = i % num_gpu

    # run script for equilibration steps : 1st --> 2nd
    if (mpi_version):
       if (set_fix):
          equil_script = '''#!/bin/sh

mpirun -np %d pmemd.MPI -O -i equil_1.in -o equil_1.out -p ../%s_solv.top -c ../%s_min2.rst -r %s_equil_1.rst -x %s_equil_1.crd -ref ../%s_min2.rst

mpirun -np %d pmemd.MPI -O -i equil_2.in -o equil_2.out -p ../%s_solv.top -c %s_equil_1.rst -r %s_equil_2.rst -x %s_equil_2.crd -ref ../%s_min2.rst

''' \
          % (num_cpu, tag, tag, tag, tag, tag, num_cpu, tag, tag, tag, tag, tag)
          if iamd != 0 :
             amd_script = 'mpirun -np %d pmemd.MPI -O -i equil_3.in -o equil_3.out -p ../%s_solv.top -c %s_equil_2.rst -r %s_equil_3.rst -x %s_equil_3.crd -ref ../%s_min2.rst'%(num_cpu, tag, tag, tag, tag, tag)
             equil_script = equil_script + amd_script
       else:
          equil_script = '''#!/bin/sh

mpirun -np %d pmemd.MPI -O -i equil_1.in -o equil_1.out -p ../%s_solv.top -c ../%s_min2.rst -r %s_equil_1.rst -x %s_equil_1.crd -ref ../%s_min2.rst

mpirun -np %d pmemd.MPI -O -i equil_2.in -o equil_2.out -p ../%s_solv.top -c %s_equil_1.rst -r %s_equil_2.rst -x %s_equil_2.crd 
''' \
          % (num_cpu, tag, tag, tag, tag, tag, num_cpu, tag, tag, tag, tag, num_cpu, tag, tag, tag, tag) 
          if iamd != 0 :
             amd_script = 'mpirun -np %d pmemd.MPI -O -i equil_3.in -o equil_3.out -p ../%s_solv.top -c %s_equil_2.rst -r %s_equil_3.rst -x %s_equil_3.crd'%(num_cpu, tag, tag, tag, tag)
             equil_script = equil_script + amd_script
    else:   # ! default : cuda version
       if (set_fix):
          equil_script = '''#!/bin/sh

export CUDA_VISIBLE_DEVICES='%d'

pmemd.cuda -O -i equil_1.in -o equil_1.out -p ../%s_solv.top -c ../%s_min2.rst -r %s_equil_1.rst -x %s_equil_1.crd -ref ../%s_min2.rst

pmemd.cuda -O -i equil_2.in -o equil_2.out -p ../%s_solv.top -c %s_equil_1.rst -r %s_equil_2.rst -x %s_equil_2.crd -ref ../%s_min2.rst
''' \
          % (gpu_dev_id, tag, tag, tag, tag, tag, tag, tag, tag, tag, tag)
          if iamd != 0 :
             amd_script = 'pmemd.cuda -O -i equil_3.in -o equil_3.out -p ../%s_solv.top -c %s_equil_2.rst -r %s_equil_3.rst -x %s_equil_3.crd -ref ../%s_min2.rst'%(tag, tag, tag, tag, tag)
             equil_script = equil_script + amd_script
       else:
          equil_script = '''#!/bin/sh

export CUDA_VISIBLE_DEVICES='%d'

pmemd.cuda -O -i equil_1.in -o equil_1.out -p ../%s_solv.top -c ../%s_min2.rst -r %s_equil_1.rst -x %s_equil_1.crd -ref ../%s_min2.rst

pmemd.cuda -O -i equil_2.in -o equil_2.out -p ../%s_solv.top -c %s_equil_1.rst -r %s_equil_2.rst -x %s_equil_2.crd
''' \
          % (gpu_dev_id, tag, tag, tag, tag, tag, tag, tag, tag, tag,)
          if iamd != 0:
             amd_script = 'pmemd.cuda -O -i equil_3.in -o equil_3.out -p ../%s_solv.top -c %s_equil_2.rst -r %s_equil_3.rst -x %s_equil_3.crd'%(tag, tag, tag, tag)
             equil_script = equil_script + amd_script
    run_equil = open('run_equil.sh','w')
    #print >> run_equil, equil_script
    #print (equil_script, end='', file=run_equil)
    run_equil.write(equil_script)
    run_equil.close()
    st = os.stat('run_equil.sh')
    os.chmod('run_equil.sh', st.st_mode | stat.S_IEXEC)
    subprocess.call('./run_equil.sh',shell=True)

    fequil1 = check_equil('equil_1.out')
    fequil2 = check_equil('equil_2.out')

    if fequil2 == 0 :
       equil_corr('equil_2',tag)
       if os.path.exists('%s_equil_2.rst'%(tag)):
          os.remove('%s_equil_2.rst'%(tag))
       if os.path.exists('%s_equil_2.out'%(tag)):
          os.remove('%s_equil_2.out'%(tag))
       if os.path.exists('%s_equil_2.crd'%(tag)):
          os.remove('%s_equil_2.crd'%(tag))
       if (set_fix):
          tequil_script = '''#!/bin/sh

export CUDA_VISIBLE_DEVICES='%d'

pmemd.cuda -O -i equil_2_t.in -o equil_2.out -p ../%s_solv.top -c %s_equil_2.rst7 -r %s_equil_2.rst -x %s_equil_2.crd -ref ../%s_min2.rst

''' \
          % (fixed_gpu_id, tag, tag, tag, tag,  tag)
          if iamd != 0 :
             amd_script = 'pmemd.cuda -O -i equil_3.in -o equil_3.out -p ../%s_solv.top -c %s_equil_2.rst -r %s_equil_3.rst -x %s_equil_3.crd -ref ../%s_min2.rst'%(tag, tag, tag, tag, tag)
             tequil_script = tequil_script + amd_script   
       else:
          tequil_script = '''#!/bin/sh

export CUDA_VISIBLE_DEVICES='%d'

pmemd.cuda -O -i equil_2_t.in -o equil_2.out -p ../%s_solv.top -c %s_equil_2.rst7 -r %s_equil_2.rst -x %s_equil_2.crd

''' \
          % (fixed_gpu_id, tag, tag, tag, tag)
          if iamd != 0 :
             amd_script = 'pmemd.cuda -O -i equil_3.in -o equil_3.out -p ../%s_solv.top -c %s_equil_2.rst -r %s_equil_3.rst -x %s_equil_3.crd'%(tag, tag, tag, tag)
             tequil_script = tequil_script + amd_script
       run_equil_t = open('run_equil_t.sh','w')
       run_equil_t.write(tequil_script)
       st = os.stat('run_equil_t.sh')
       os.chmod('run_equil_t.sh', st.st_mode | stat.S_IEXEC)
       subprocess.call(['./run_equil_t.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
       fequil2 = check_equil('equil_2.out')
       fequil3 = check_equil('equil_3.out')
    else:
       fequil3 = check_equil('equil_3.out')
    ###############################################
    ##---- production dir
    if fequil1 > 0 and fequil2 > 0 and fequil3 > 0 :
       os.chdir('../%s' % (prod_dir))
       if (set_fix):
          shutil.copy('%s/simulation_in/tmp_md_fix.in' % (common), './md.in')
          set_sim_options('./md.in', process='PROD', str_res_id=lstart, end_res_id=lres) 
       else:
          shutil.copy('%s/simulation_in/tmp_md.in' % (common), './md.in')
          print(os.getcwd())
          set_sim_options('./md.in', process='PROD')
    
       if (mpi_version):
          shutil.copy('%s/md_local_script/prep_production_run_MPI.py' % (bin_dir), './prep_production_run_MPI.py')
          os.system('./prep_production_run_MPI.py -np %d -t %d -f %d -n %s -p ../%s_solv.top -r %s_equil_2.rst' \
                 % (num_cpu, \
                    # num of cpu for MPI
                    int(prod_total_time), \
                    # total production time (overall time, ns)
                    int(prod_freq_time * 1000), \
                    # one production time ( ps )
                    '%s_%dK_%s' % (tag, int(target_temp), traj_id),\
                    # name of this simulation trajectory
                    tag, tag))     # name of topology
       else:
          shutil.copy('%s/md_local_script/prep_production_run_CUDA_v1.py' % (bin_dir), './prep_production_run_CUDA.py')
          if (set_fix):
             if iamd != 0:
                os.system('./prep_production_run_CUDA.py -i %d -t %d -f %d -n %s -p ../%s_solv.top -r %s_equil_3.rst -fix'%(gpu_dev_id,int(prod_total_time),int(prod_freq_time * 1000),'%s_%dK_%s' % (tag, int(target_temp), traj_id),tag, tag))   
             else:
                os.system('./prep_production_run_CUDA.py -i %d -t %d -f %d -n %s -p ../%s_solv.top -r %s_equil_2.rst -fix' \
                       % (gpu_dev_id, \
                           # GPU card id
                           int(prod_total_time), \
                           # total production time (overall time, ns)
                           int(prod_freq_time * 1000), \
                           # one production time ( ps )
                           '%s_%dK_%s' % (tag, int(target_temp), traj_id),\
                           # name of this simulation trajectory
                           tag, tag))
          else:
             if iamd != 0:
                os.system('./prep_production_run_CUDA.py -i %d -t %d -f %d -n %s -p ../%s_solv.top -r %s_equil_3.rst'%(gpu_dev_id,int(prod_total_time),int(prod_freq_time * 1000),'%s_%dK_%s' % (tag, int(target_temp), traj_id),tag, tag))
             else:
                os.system('./prep_production_run_CUDA.py -i %d -t %d -f %d -n %s -p ../%s_solv.top -r %s_equil_2.rst' \
                       % (gpu_dev_id, \
                          # GPU card id
                          int(prod_total_time), \
                          # total production time (overall time, ns)
                          int(prod_freq_time * 1000), \
                          # one production time ( ps )
                          '%s_%dK_%s' % (tag, int(target_temp), traj_id),\
                          # name of this simulation trajectory
                          tag, tag))     # name of topology
       subprocess.call('./run_production.sh',shell=True)
    
       os.chdir('../../')
