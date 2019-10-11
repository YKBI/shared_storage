import os,sys
import pandas as pd
import numpy as np
import subprocess as sp
import shutil
#argument
tag     = sys.argv[1] # out directory tag
pdbna   = sys.argv[2] # pdb name
seq     = sys.argv[3] # sequence
HLAcl   = sys.argv[4] # HLA class example == 'A' 'B' 'C'
HLAty   = sys.argv[5] # HLA type example == '0201' '0206'
num_cpu = sys.argv[6] # number of cpu for docking
nst     = sys.argv[7] # number of conformation (number of docking out pdbs) default is 1000

#directories
ROOT        = '/awork06-1/YKLee/'
FLEXPEP_BIN = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/source/bin'
ROSETTA_DB  = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/main/database'
GEAR        = '/awork08/93_hong/NGS_ARS/Gear'
REC_dir     = ROOT + 'receptor/' + HLAcl + '/' + HLAty + '/'
PEP_dir     = ROOT + 'peptide/' + HLAcl + '/' + HLAty + '/'
WD          = seq + '_' + pdbna + '_' + HLAcl + HLAty
wdir = os.getcwd()
#make directory structure
dir_list = [REC_dir,PEP_dir,WD,WD+'/PDB_1ST',WD+'/PDB_2ND',WD+'/DOCK_RES']
for dd in dir_list:
    try:
        if not os.path.exists(dd):os.makedirs(dd,0777)
    except OSError:pass

#clean_pdb.py by Rosetta
def clean(x):
    cmd1 = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py ' + x + '.pdb A'
    cmd2 = '/awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py ' + x + '.pdb B'
    sp.call(cmd1,shell=True)
    sp.call(cmd2,shell=True)
    shutil.move(x + '_A.pdb',REC_dir + x + '_rec.pdb')
    shutil.move(x + '_B.pdb',PEP_dir + x + '_pep.pdb')


#make up peptide
def ending_pepoc(x,y):
    mutation = 'mutate.pl -seq 1:'  + y + ' ' + PEP_dir + x + '_pep.pdb > ' + x + '_' + y + '_temp.pdb'
    sp.call(mutation,shell=True)
    temp_list = [PEP_dir + x + '_pep.pdb',x + '_' + y + '_temp.pdb']
    for ppdb in temp_list:
        with open(ppdb,'r') as F:
            Lines = []
            if 'temp' in ppdb:
                for line in F.readlines():
                    if line[21:23].strip() != 'B' :sLines =line
                    elif line[21:23].strip() == 'B' or line[21:23].strip() == '':sLines ==lin[:56] + '1' + line[57:]
                Lines.append(sLines)
            else:
                for line in F.readlines():
                    if line[21:23].strip() != 'B' or line[21:23].strip() == '': sLines = line
                    elif line[21:23].strip() =='B': sLines = line[:56] + '1' + line[57:]
                Lines.append(sLines)
        with open(ppdb + '_new.pdb','w') as N:
            for LineOut in Lines:N.write(LineOut)
            N.write('END')
        shutil.move(ppdb + '_new.pdb',ppdb)

#make up receptor
def helpa(x):
    temp_list = []
    #def n2
    with open(x + '.pdb.st','r') as S:
        for line in S.readlines():
            if line.startswith('LOC') and line.find('AlphaHelix') >0:
                temp_list.append(int(line.split()[6]))
                n2 = int(max(temp_list)) + 3
    #cutting receptor
	REC_Lines = []
    with open(REC_dir + x + '_rec.pdb','r') as F:
        for line in F.readlines():
            if line.startswith('ATOM'):
                if 1 <= int(line[22:26].strip()) <= n2 : REC_Lines.append(line)
                else : continue
            else : REC_Lines.append('TER\n')
    with open(x + 'new.pdb','w') as N:
        for line in REC_Lines: N.write(line)
    shutil.move(x + 'new.pdb', REC_dir + x +'_rec.pdb')

def gen_inp(x,y):
    cmd = 'cat ' + REC_dir + x + '_rec.pdb ' + x + '_' + y + '_temp.pdb > ' + x + '_' + y + '_inp.pdb'
    sp.call(cmd,shell=True)
#generate native pdb (no mutation)
def gen_native(xx,yy,zz):
    HLAclaty = xx[-1] + yy
    OLines = []
    with open(REC_dir + zz + '_rec.pdb','r') as rpdb:
        for line in rpdb.readlines() : OLines.append(line)
        OLines.append('TER\n')
    with open(PEP_dir + zz + '_pep.pdb','r') as ppdb:
        for line in ppdb.readlines():
            if line.startswith('ATOM') and line.find('OXT') < 0:
                sLine = line[:21] + 'B ' + line[23:]
                OLines.append(sLine)
        OLines.append('END')
    try:
        if not os.path.exists(ROOT + HLAclaty + '_native'): os.mkdir(ROOT + HLAclaty + '_native')
    except OSError: pass
    with open(ROOT + HLAclaty + '_native/' + zz + '_native.pdb','w') as npdb:
        for line in OLines:npdb.write(line)

def sheba_run(x,y):
    sheba1 = GEAR + '/sheba_01 -x ' + x + '.pdb ' + y + '.pdb'
    sheba2 = GEAR + '/sheba_01 -t ' + y + '.trf ' + y + '.pdb'
    sp.call(sheba1,shell=True)
    sp.call(sheba2,shell=True)
    shutil.move(y + '.pdb.pdb',y + '_tr.pdb')

def ext_pep(x,y):
    OLines = []
    with open(x + '.pdb','r') as F:
        for line in F.readlines():
            if line.startswith('ATOM') and line[21] == y:OLines.append(line)
        OLines.append('END\n')
    with open(x + '_pep.pdb','w') as F:
        for line in OLinese:F.write(line)

def docking(x,y,z,n):
    HLAclaty = HLAcl[-1] + HLAty
    prepack_Lines = ['-s ' + x + '_' + y + '_inp.pdb',\
                     '-out:path:all ' +WD,\
                     '-out:prefix PPK_',\
                     '-out:file:scorefile score_ppk_' + x + '_' + y + '.sc'\
                     ,'-ex1',\
                     '-ex2aro',\
                     '-use_input_sc',\
                     '-flexpep_prepack',\
                     '-nstruct 1']
    flag1_Lines = ['-s ' + WD + '/PPK_' + x + '_' + y + '_inp_0001.pdb',\
                  '-native ' + ROOT + HLAclaty + '_native/' + x + '_native.pdb',\
                  '-out:path:all ' + WD + '/PDB_1ST',\
                  '-out:file:scorefile score_1st_' + x + '_' + y + '.sc',\
                  '-nstruct ' + z,\
                  '-flexPepDocking:lowres_preoptimize',\
                  '-flexPepDocking:pep_refine',\
                  '-flexPepDocking:flexpep_score_only',\
                  '-ex1',\
                  '-ex2aro']
    flag2_Lines = ['-s ' + WD + '/PPK_' + x + '_' + y + '_inp_0001.pdb', \
                   '-native ' + ROOT + HLAclaty + '_native/' + x + '_native.pdb', \
                   '-out:path:all ' + WD + '/PDB_2ND', \
                   '-out:file:scorefile score_1st_' + x + '_' + y + '.sc', \
                   '-nstruct ' + z, \
                   '-flexPepDocking:lowres_preoptimize', \
                   '-flexPepDocking:pep_refine', \
                   '-flexPepDocking:flexpep_score_only', \
                   '-ex1', \
                   '-ex2aro']
    with open(WD + "/prepack_flags",'w') as prepack:
        for line in '\n'.join(prepack_Lines):prepack.write(line)
    cmd = FLEXPEP_BIN + '/FlexPepDocking.linuxgccrelease -database ' + ROSETTA_DB + ' @' + WD + '/prepack_flags > ' + WD + '/prepack.log'
    sp.call(cmd,shell=True)
    with open(WD + '/run_flags1' , 'w') as f1:
        for line in '\n'.join(flag1_Lines):f1.write(line)
    with open(WD + '/run_flags2' , 'w') as f2:
        for line in '\n'.join(flag2_Lines):f2.write(line)

    cmd1 = 'mpirun -np ' + n + ' ' + FLEXPEP_BIN + '/FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + WD + '/run_flags1 > ' + WD + '/run1.log'
    cmd2 = 'mpirun -np ' + n + ' ' + FLEXPEP_BIN + '/FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + WD + '/run_flags2 > ' + WD + '/run2.log'
    sp.call(cmd1,shell=True)
    sp.call(cmd2,shell=True)

def make_tsv(x,y):
    temp_list = []
    with open(WD + '/PDB_1ST/score_1st_' + x + '_' + y + '.sc','r') as F:
        for line in F.readlines():
            if line.startswith('SCORE:') :
                line_1 = ' '.join(line.split())
                col = line_1.split(' ')
                if line_1.find('total_score'): temp_list.append(col[68] + '\t' + col[45] + '\t' + col[1])
                else : temp_list.append('PDB_1ST/' + col[68] + '.pdb\t' + col[45] + '\t' + col[1])
    with open(WD + '/PDB_2ND/score_2nd_' + x + '_' + y + '.sc','r') as F:
        for line in F.readlines():
            if line.startswith('SCORE:') and line.find('total_score') < 0:
                line_1 = " ".join(line.split())
                col = line_1.split(' ')
                temp_list.append('PDB_2ND/' + col[62] + '.pdb\t' + col[45] + '\t' + col[1] )
    with open(WD + '/total_score.tsv','w') as TS:
        for line in '\n'.join(temp_list): TS.write(line)

    df = pd.read_csv(WD + '/total_score.tsv' ,sep = '\t').sort_values(['total_score']).head(10)
    df.to_csv(WD + '/dock_selected.tsv',sep='\t',index=False)
    for line in df['description']:
        if line.split('/')[0].find('1ST'):shutil.copy(WD + '/' + line, WD + '/DOCK_RES/1ST_' + line.split('/')[1])
        elif line.split('/')[0].find('2ND'):shutil.copy(WD + '/' + line, WD + '/DOCK_RES/2ND_' + line.split('/')[1])

def cal_rms(x):
    pdbs = []
    with open(ROOT + HLAcl + '_list/' + HLAty + '.list','r') as List:
        for line in List.readlines():pdbs.append(line.strip())
    for ipdb in pdbs:
        if not os.path.exists(ROOT + HLAcl + '_native/' + ipdb + '_native.pdb'):
            gen_native(HLAcl,HLAty,ipdb)
        shutil.copy(ROOT + HLAcl + '_native/' + ipdb + '_native.pdb',WD + '/DOCK_RES')
    os.chdir(WD + '/DOCK_RES')
    dock_list = glob.glob('*.pdb')
    for ipdb in pdbs :
        if ipdb == x:
            if not os.path.exists(ipdb + '_native_pep.pdb'):ext_pep(x + '_native','B')
            for dock in dock_list:
                if not os.path.exists(dock.split('.')[0] + '_pep.pdb'):ext_pep(dock.split('.')[0],'B')
                rmsd = GEAR + '/rmsd_total ' + ipdb +'_native_pep.pdb ' + dock.split('.')[0] + '_pep.pdb bb >> rmsd_total.txt'
                sp.call(rmsd,shell=True)
        else:
            sheba_run(x + '_native',ipdb + '_native')
            ext_pep(ipdb + '_native_tr','B')
            for dock in dock_list:
                if not os.path.exists(dock.split('.')[0] + '_pep.pdb'):
                    ext_pep(dock.split('.')[0],'B')
                rmsd = GEAR + '/rmsd_total ' + ipdb + '_native_tr_pep.pdb ' + dock.split('.')[0] + '_pep.pdb bb >> rmsd_total.txt'
                sp.call(rmsd,shell=True)

    df2 = pd.read_csv('rmsd_total.txt',sep='\t',names =('DOCK_STR','HLA_REF','PEP_RMSD'))
    df2.sort_values(['PEP_RMSD']).to_csv('rmsd_total_sort.txt',sep='\t',index=False)

#doing clean ~ docking
clean(pdbna)
ending_pepoc(pdbna,seq)
helpa(pdbna)
if not os.path.exists(ROOT + HLAcl + '_native/' + pdbna + '_native.pdb'):
    gen_native(HLAcl,HLAty,pdbna)
gen_inp(pdbna,seq)
docking(pdbna,seq,nst,num_cpu)
#calculate rmsd
make_tsv(pdbna,seq)
os.chdir(wdir)
cal_rms(pdbna)
os.chdir(wdir)
