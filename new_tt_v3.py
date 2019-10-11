import os,sys,glob,shutil,string
import subprocess as sp
import pandas as pd
GEAR = '/awork06-1/neoscan_gear'
inp_dir = sys.argv[1]
tag = sys.argv[2]
out_arg = sys.argv[3]


def remove_duplicates(x):
    my_set = set()
    res = []
    for e in x:
        if e not in my_set:
            res.append(e)
            my_set.add(e)
    return res

def help():
    print "Usage: python %s [inp_dir] [PDB] [output folder]\n" %sys.argv[0]
    return

if len(sys.argv) == 1:
    help()
else:
    try:
        if not os.path.exists(out_arg + '/orig_pdb'):
            os.makedirs(out_arg + '/orig_pdb',0777)
    except OSError: pass
    try:
        if not os.path.exists(out_arg + '/score'):
            os.makedirs(out_arg + '/score',0777)
    except OSError: pass

    os.chdir(inp_dir + '/traj_1/pdb_from_prod')
    sfpdbs = glob.glob('*.pdb')
    fpdbs = []
    for pdb in sfpdbs:
        if '_tr' in pdb: pass
        else: fpdbs.append(pdb)
    for fpdb in sorted(fpdbs):
        ser = int(fpdb.split('.')[2])
        nfpdb = tag + '.' + str(ser)
        try:
            if not os.path.exists(out_arg + '/' + nfpdb):
                os.makedirs(out_arg + '/' + nfpdb,0777)
        except OSError: pass
        idx = 0
        with open(out_arg + '/orig_pdb/' + nfpdb + '.pdb','w') as f:
            with open(fpdb,'r') as f1:
                for line in f1.readlines():
                    if line.startswith('ATOM') > 0 and line.find('OXT') < 0:
                        if idx == 0 and line[21] != 'A':
                            TEXT = line[:21] + 'A' + line[22:]
                            f.write(TEXT)
                        elif idx > 0 and line[21] != 'B':
                            TEXT = line[:21] + 'B' + line[22:]
                            f.write(TEXT)
                        else: f.write(line)
                    elif line.startswith('TER') > 0:
                        if idx == 0: f.write(line)
                        idx += 1
            f.write('END')

    	cmd = GEAR + '/reduce -Trim ' + out_arg + '/orig_pdb/' + nfpdb + '.pdb > ' + out_arg + '/orig_pdb/' + nfpdb + '_red.pdb'
    	sp.call(cmd,shell=True)

    	with open(out_arg + '/' + nfpdb + '/' + nfpdb + '_rec1.pdb','w') as rf:
        	with open(out_arg + '/orig_pdb/' + nfpdb +'_red.pdb','r') as of:
            		for line in of.readlines():
                		if line.startswith('ATOM') > 0 and line[21] =='A': rf.write(line)
        		rf.write('TER\n')

    	with open(out_arg + '/' + nfpdb + '/' + nfpdb + '_pep.pdb','w') as pf:
        	with open(out_arg + '/orig_pdb/' + nfpdb + '_red.pdb','r') as of:
            		for line in of.readlines():
                		if line.startswith('ATOM') > 0 and line[21] == 'B': pf.write(line)
        		pf.write('END')
    tf = pd.read_csv('ext_bb_rmsd.dat',sep='\t',header=None)
    tf.columns = ['PDB','rmsBB']
    tf.to_csv(out_arg + '/score/' + pdb + '_total_score.tsv',header=True,index=False,sep='\t')
