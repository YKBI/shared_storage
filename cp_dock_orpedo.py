import sys
import os
import shutil
import subprocess as sp

inp_dir = sys.argv[1] #input directory
num_cpu = sys.argv[2] #N. of CPU
#pdb = sys.argv[3]
pre_dir = os.getcwd()
seq = inp_dir.split('_')[0]
pdb = inp_dir.split('_')[1]
HLAclass = inp_dir.split('_')[2]
HLAtype = inp_dir.split('_')[3]
HLA = '_'.join(inp_dir.split('_')[1:])
# over 0.6 r-square score pdbs
'''
pdb_dic = {'HLA-A_0201':['1oga','2gtw','3d25','3mrb','5eu5'],\
          'HLA-A_0206':['3mre','3mrk','3mrg','3mrd','3gso'],\
          'HLA-A_1101':['1q94','1x7q','4n8v','2uwe','2av7'],\
          'HLA-A_2402':['2bck','2jcc','3i6l','4nqx','3bo8'],\
          'HLA-A_3303':['4hx1','4hwz','2xpg','3rl1','4i48'],\
          'HLA-B_1501':['1xr8','5txs','3c9n','1a9e','1xr9'],\
          'HLA-B_3501':['2cik','6bj8','3lko','3lkp','3lkq'],\
          'HLA-B_4403':['1n2r','1sys','3kpn','3kpo','3l3i'],\
          'HLA-B_5101':['1a1m','1a1o','3c9n','3x14','1e27'],\
          'HLA-B_5401':['3bvn','3bxn','4u1j','4u1m','1mi5'],\
           'HLA-C_0102':['5w69','5w6a'],\
	   'HLA-C_0303':['1qqd','5xs3'],\
           'HLA-C_0304':['1efx','1im9'],\
           'HLA-C_0702':['5vge','5w67'],\
           'HLA-C_0801':['5vgd','4nt6']}
'''
#Change pdb's occupancy
def pepoc(x):
	with open(x,'r') as F:
		Lines =[]
		for line in F.readlines():
			if line[21:23].strip() != 'B' or line[21:23].strip() =='': sLines=line
			elif line[21:23].strip() =='B': sLines=line[:56]+'1'+line[57:]
			Lines.append(sLines)
	with open('%s_new.pdb'%x.split('.')[0],'w') as N:
		for LineOut in Lines:
			N.write(LineOut)
	shutil.move('%s_new.pdb'%x.split('.')[0],x)
#Operation docking scripts
def dock(pdb):
    docking = 'python /awork06-1/YKLee/run_flexpep_dock_v1.py %s %s %s %s %s'%(seq,pdb,HLAclass,HLAtype,num_cpu)
    sp.call(docking,shell=True)

shutil.copytree(pre_dir+'/'+inp_dir,'/home/user1/'+inp_dir)
os.chdir('/home/user1/'+inp_dir)
'''
for pdb in pdb_dic[HLA]:
    pepoc(pdb+'_'+seq+'_inp.pdb')
    dock(pdb)
'''
pepoc(pdb+'_'+seq+'_inp.pdb')
dock(pdb)
os.chdir('../')
zip = 'zip -r ' + inp_dir + '.zip ' + inp_dir
sp.call(zip,shell=True)
'''
try:
    os.mkdir('/awork06-1/YKLee/orpedo_out/')

except:
    shutil.copytree(inp_dir,'/awork06-1/YKLee/orpedo_out/'+inp_dir)
shutil.copytree(inp_dir,'/awork06-1/YKLee/orpedo_out/'+inp_dir)
'''
