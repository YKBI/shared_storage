import sys,os
import subprocess as sp
import glob
inp_dir = sys.argv[1]
nst = sys.argv[2]
ncpu = sys.argv[3]
out_arg = sys.argv[4]
ref = '_'.join(inp_dir.split('_')[:2])
os.chdir(inp_dir)
lig_dic = {}
def make_dic(x):
    pdb_list = sorted(glob.glob('%s_*'%x))
    for pdb in pdb_list:
        lig_dic['_'.join(pdb.split('_')[:5])] = pdb.split('_')[-2]
def docking(x,y,z,n):
    docking = 'time python /awork06-1/YKLee/run_lig_dock_v3.py %s %s %s %s '%(x,y,z,n)
    sp.call(docking,shell=True)
try:
    cpc = 'cp -r ' + inp_dir + ' /home/user1/' + inp_dir
    sp.call(cpc,shell=True)
except: pass
os.chdir('/home/user1/'+inp_dir)
make_dic(ref)
for lig in lig_dic:
    docking(lig.split('.')[0],lig_dic[lig],nst,ncpu)
os.chdir('../')
tar = 'tar -zcvf ' + inp_dir + '.tar.gz ' + inp_dir
sp.call(tar,shell=True)
try:
    if not os.path.exists('/awork05-1/YKLee/kinase_dock/' + out_arg):
        os.makedirs('/awork05-1/YKLee/kinase_dock/' + out_arg)
except OSError: pass
try:
    cpcc  = 'cp -r ' + inp_dir + '/awork05-1/YKLee/kinase_dock/' + out_arg + '/' + inp_dir
    sp.call(cpcc,shell=True)
except : pass