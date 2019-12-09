import sys,os,subprocess,glob

inp_dir = sys.argv[1]
os.chdir(inp_dir+'/dock_res/')
GEAR = '/awork06-1/neoscan_gear'
def enva_run(i):
	ff=i.split('.')[0]
	enva1 = GEAR+'/enva_rec3 -c ' + ff + '.pdb > '+ff + '.out'
#	enva2 = GEAR+'/enva_rec3 -b ' + ff + '.pdb > '+ff + '_bb.out'
#	enva3 = GAER+'/enva_rec3 -a ' + ff + '.pdb > '+ff + '_aa.out'
	
	return enva1
def run_ep(p):
	proc = subprocess.Popen([GEAR+'/enva_rec3 -c',p])
	return proc
procs = []
for i in sorted(glob.glob('*.pdb')):
	proc = run_ep(i)
	procs.append(proc)
