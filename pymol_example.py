import subprocess as sp
import sys,glob,os

pdbna = sys.argv[1]
pdb_list = []
#make pdb list
with open('list','r') as F:
    for pdb in F.readlines():
        pdb_list.append(pdb.strip())

def ext(x):
    pml_list = []
    pdbs =  sorted(glob.glob(x.strip()+'_trans/'+ x.strip() + '_*.pdb'))
    for pdb in pdbs:
        pml_list.append('load ' + pdb.strip())
    pp_list=['bg_color white','color green,*','color red,het','show cartoon,*','show sticks,het','orient','hide line','rotate [0,0,1],90','set seq_view,1','set seq_view_label_mode,3']
    with open(x.strip()+'.pml','w') as W:
        for line in pml_list:
            W.write(line + '\n')
        for lin in pp_list:
            W.write(lin + '\n')

    cmd = 'pymol '+ x.strip() +'.pml'
    sp.call(cmd,shell=True)
ext(pdbna)
