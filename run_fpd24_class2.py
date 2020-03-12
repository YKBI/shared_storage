import os,sys,psutil,glob
import subprocess
import pandas as pd
seq = sys.argv[1]
pdb = sys.argv[2]
HLAclass = sys.argv[3]
HLAtype = sys.argv[4]
num_cpu = str(psutil.cpu_count()-1)
nout = sys.argv[5]
switch = sys.argv[6]
HLAclaty = HLAclass[-1] + HLAtype
RES = seq + '_' + pdb + '_' + HLAclass + HLAtype
ROOT = '/home/user1/'
os.environ['rosetta'] = '/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin/'
rosetta = os.environ['rosetta']
os.environ['PATH'] += ':' + rosetta + ':' + os.environ['PATH']
#FLEXPEP_BIN = '/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin'
ROSETTA_DB = '/lwork01/rosetta_src_2019.40.60963_bundle/main/database'
GEAR = '/home/user1/neoscan_gear'
PEP_cwd = '/home/user1/minimize/'
REC_cwd = '/home/user1/minimize/'
Olines = []
def mutate(x,y):
	with open('%s_%s_tmp.pdb'%(x,y),'r') as F:
		for line in F.readlines():
			if not line.startswith('TER\n' or 'END\n'):
				Olines.append(line[:21] + 'B ' + line[23:56] + '1' + line[57:])
			else:
				Olines.append(line)
		del Olines[-1]
	with open('%s_%s_t1.pdb'%(x,y),'w') as W:
		W.write('TER\n')
		for line in Olines:
			W.write(line)
	os.system('cat %s_rec.pdb %s_%s_t1.pdb > %s_%s_inp.pdb'%(x,x,y,x,y))
				
def gen_native(aa,bb,cc):
    HLAclaty = aa[-1] + bb
    OLines = []
    #with open(REC_cwd + RES + '/' + cc + '_rec.pdb','r') as rpdb:
    with open(REC_cwd + '/' + cc + '_rev_A.pdb','r') as rpdb:
        for line in rpdb.readlines(): OLines.append(line)
        OLines.append('TER\n')
    #with open(PEP_cwd + RES + '/'+ cc + '_pep.pdb','r') as ppdb:
    with open(PEP_cwd + '/' + cc + '_rev_B.pdb','r') as ppdb:
        for line in ppdb.readlines():
            if line.startswith('ATOM') and line.find('OXT') <0 :
                sline = line[:21] + 'B ' + line[23:]
                OLines.append(sline)
        OLines.append('END\n')
    try:
        if not os.path.exists(ROOT + 'Native/' + HLAclaty + '_native'):
            os.makedirs(ROOT + 'Native/' + HLAclaty + '_native')
    except OSError:
        pass
    with open(ROOT + 'Native/' + HLAclaty + '_native/' + cc + '_native.pdb','w') as npdb:
        for line in OLines: npdb.write(line)

def sheba_run(ref,org):
    sheba1 = GEAR + '/sheba_01 -x ' + ref + '.pdb ' + org + '.pdb'
    sheba2 = GEAR + '/sheba_01 -t ' + org + '.trf ' + org + '.pdb'
    subprocess.call(sheba1,shell=True)
    subprocess.call(sheba2,shell=True)
    shutil.move(org + '.pdb.pdb',org + '_tr.pdb')

def ext_pep(pdbx,ch):
    OLines = []
    with open(pdbx +'.pdb','r') as F:
        for line in F.readlines():
            if line.startswith('ATOM') and line[21] == ch: OLines.append(line)
        OLines.append('END\n')
    with open(pdbx+'_pep.pdb','w') as F:
        for line in OLines: F.write(line)

wdir = os.getcwd()
if int(switch) == 1:
	mutate(pdb,seq)
else:
	pass

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

if not os.path.exists(ROOT + '/Native/' + HLAclaty + '_native/' + pdb + '_native.pdb'):
	gen_native(HLAclass,HLAtype,pdb)

with open(RES + "/prepack_flags",'a') as pf :
    pf.write('-s ' + pdb + '_' + seq + '_inp.pdb\n')
    pf.write('-out:path:all ' + RES + '\n')
    pf.write('-out:prefix PPK_\n')
    #	pf.write('-out:suffix test\n')
    pf.write('-out:file:scorefile score_ppk_' + pdb + '_' + seq + '.sc\n')
    pf.write('-ex1\n')
    pf.write('-ex2aro\n')
    pf.write('-use_input_sc\n')
    pf.write('-flexpep_prepack\n')
    pf.write('-nstruct 1\n')

os.system('FlexPepDocking.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/prepack_flags > ' + RES + '/prepack.log')

with open(RES + "/run_flags1",'a') as rf :
    rf.write('-s ' + RES + '/PPK_' + pdb + '_' + seq + '_inp_0001.pdb\n')
    rf.write('-native ' + ROOT + 'Native/' + HLAclaty + '_native/' + pdb + '_native.pdb\n')
    rf.write('-out:path:all ' + RES + '/PDB_1ST\n')
    rf.write('-out:file:scorefile score_1st_' + pdb + '_' + seq + '.sc\n')
    rf.write('-nstruct %s\n'%nout)
    rf.write('-flexPepDocking:lowres_preoptimize\n')
    rf.write('-flexPepDocking:pep_refine\n')
    rf.write('-flexPepDocking:flexpep_score_only\n')
    rf.write('-ex1\n')
    rf.write('-ex2aro\n')

with open(RES + "/run_flags2",'a') as rf1 :
    rf1.write('-s ' + RES + '/PPK_' + pdb + '_' + seq + '_inp_0001.pdb\n')
    rf1.write('-native ' + ROOT + 'Native/' + HLAclaty + '_native/' + pdb + '_native.pdb\n')
    rf1.write('-out:path:all ' + RES + '/PDB_2ND\n')
    rf1.write('-out:file:scorefile score_2nd_' + pdb + '_' + seq + '.sc\n')
    rf1.write('-nstruct %s\n'%nout)
    rf1.write('-flexPepDocking:pep_refine\n')
    rf1.write('-flexPepDocking:flexpep_score_only\n')
    rf1.write('-ex1\n')
    rf1.write('-ex2aro\n')

os.system('mpirun -np ' + num_cpu + ' FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/run_flags1 > ' + RES + '/run1.log')
os.system('mpirun -np ' + num_cpu + ' FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/run_flags2 > ' + RES + '/run2.log')

with open(RES + '/total_score.tsv','a') as fx :
    with open(RES + '/PDB_1ST/score_1st_' + pdb + '_' + seq + '.sc','r') as f :
        lines = f.readlines()
        for line in lines:
            if line.startswith('SCORE:') > 0 :
                line = " ".join(line.split())
                col = line.split(' ')
                if line.find('total_score') > 0 :
                    if len(seq) == 8:
                        fx.write('%s\t%s\t%s\n'%(col[67],col[44],col[1]))#(col[68],col[45],col[1]))
		    elif len(seq) == 9:
			fx.write('%s\t%s\t%s\n'%(col[68],col[45],col[1]))
		    elif len(seq) == 10:
			fx.write('%s\t%s\t%s\n'%(col[69],col[46],col[1]))
                else :
		    if len(seq) == 8:
                        fx.write('PDB_1ST/%s.pdb\t%s\t%s\n'%(col[67],col[44],col[1]))#(col[68],col[45],col[1]))
		    elif len(seq) == 9:
			fx.write('PDB_1ST/%s.pdb\t%s\t%s\n'%(col[68],col[45],col[1]))
		    elif len(seq) == 10:
			fx.write('PDB_1ST/%s.pdb\t%s\t%s\n'%(col[69],col[46],col[1]))
    with open(RES + '/PDB_2ND/score_2nd_' + pdb + '_' + seq + '.sc','r') as f1 :
        lines = f1.readlines()
        for line in lines :
            if line.startswith('SCORE:') > 0 and line.find('total_score') < 0 :
                line = " ".join(line.split())
                col = line.split(' ')
		if len(seq) == 8:
                    fx.write('PDB_2ND/%s.pdb\t%s\t%s\n'%(col[61],col[44],col[1]))
		elif len(seq) == 9:
		    fx.write('PDB_2ND/%s.pdb\t%s\t%s\n'%(col[62],col[45],col[1]))
		elif len(seq) == 10:
		    fx.write('PDB_2ND/%s.pdb\t%s\t%s\n'%(col[63],col[46],col[1]))



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

with open(ROOT + HLAclaty + '_' + str(len(seq)) + 'seq.list','r') as f3 :
    if len(seq) != 9:
        for line in f3.readlines() : pdbs.append(line.split()[2].lower())
    else:
	for line in f3.readlines() : pdbs.append(line.strip().lower())

for ipdb in pdbs :
    if not os.path.exists(ROOT + 'Native/' + HLAclaty + '_native/' + ipdb + '_native.pdb'):
        gen_native(HLAclass, HLAtype, ipdb)
    shutil.copy(ROOT + 'Native/' + HLAclaty + '_native/' + ipdb + '_native.pdb',RES + '/DOCK_RES')

os.chdir(RES + '/DOCK_RES')
for ipdb in pdbs :
    if pdb == ipdb :
        if not os.path.exists(ipdb + '_native_pep.pdb') :
            ext_pep(pdb + '_native','B')
        for dock in dock_list :
            if not os.path.exists(dock.split('.')[0] + '_pep.pdb'):
                ext_pep(dock.split('.')[0],'B')
            cmd4 = GEAR + '/rmsd_total ' + ipdb + '_native_pep.pdb ' + dock.split('.')[0] + '_pep.pdb bb >> rmsd_total.txt'
            subprocess.call(cmd4,shell=True)
    else :
        sheba_run(pdb + '_native',ipdb + '_native')
        ext_pep(ipdb + '_native_tr','B')
        for dock in dock_list :
            if not os.path.exists(dock.split('.')[0] + '_pep.pdb'):
                ext_pep(dock.split('.')[0],'B')
            cmd4 = GEAR + '/rmsd_total ' + ipdb + '_native_tr_pep.pdb ' + dock.split('.')[0] + '_pep.pdb bb >> rmsd_total.txt'
            subprocess.call(cmd4,shell=True)

df2 = pd.read_csv('rmsd_total.txt',sep='\t',names=('DOCK_STR','HLA-REF','PEP_RMSD'))
df2 = df2.sort_values(['PEP_RMSD'])
df2.to_csv('rmsd_total_sort.txt',sep='\t',index=False)
os.chdir(wdir)
