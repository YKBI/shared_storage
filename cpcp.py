import sys,os
import subprocess as sp
import shutil
inp_dir = sys.argv[1]
o_dir_arg = sys.argv[2]
num_cpu = sys.argv[3]
pre_dir = os.getcwd()
seq = inp_dir.split('_')[0]
HLAclass = inp_dir.split('_')[1]
HLAtype = inp_dir.split('_')[2]
HLA = '_'.join(inp_dir.split('_')[1:])
HH = ''.join(inp_dir.split('_')[1:])
o_dir = '/awork05-1/YKLee/' + o_dir_arg + '/' + seq
pdb_dic = {'HLA-A_0201':['1oga','2gtw'],\
          'HLA-A_0206':['3mre','3mrk'],\
          'HLA-A_1101':['1q94','1x7q'],\
          'HLA-A_2402':['2bck','2jcc','3i6l'],\
          'HLA-A_3303':['4hx1'],\
          'HLA-B_1501':['1xr8','5txs'],\
          'HLA-B_3501':['2cik'],\
          'HLA-B_4403':['1n2r','1sys','3kpn','3kpo'],\
          'HLA-B_5101':['1a1m','1a1o','3c9n'],\
          'HLA-B_5401':['3bvn','3bxn','4u1j','4u1m'],\
           'HLA-C_0303':['1qqd'],\
           'HLA-C_0304':['1efx','1im9'],\
           'HLA-C_0702':['5vge'],\
           'HLA-C_0801':['5vgd']}

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
	docking = 'time python /awork06-1/YKLee/run_flexpep_dock_v1.py %s %s %s %s %s'%(seq,pdb,HLAclass,HLAtype,num_cpu)
	sp.call(docking,shell=True)
def tt_v3(seq,pdb,HH,o_dir):
	cmd1 = 'time python /awork06-1/YKLee/tt_v3.py %s %s %s %s'%(seq,pdb,HH,o_dir)
	sp.call(cmd1,shell=True)
try:
	cpc =  'cp -r ' + inp_dir + ' /home/user1/' +inp_dir
	sp.call(cpc,shell=True)
except:
	pass
os.chdir('/home/user1/'+inp_dir)
for pdb in pdb_dic[HLA]:
	pepoc(pdb+'_'+seq+'_inp.pdb')
	dock(pdb)
	print(os.getcwd())
	tt_v3(seq,pdb,HH,o_dir)
os.chdir('../')
try:
	if not os.path.exists(o_dir):
		os.makedirs(o_dir)
except OSError:
	pass
#for pdb in pdb_dic[HLA]:
#	tt_v3(seq,pdb,HLA,o_dir)
#os.chdir('../')
zip1 = 'time zip -r ' + inp_dir + '.zip ' + inp_dir
sp.call(zip1,shell=True)
cpcp = 'time cp ' + inp_dir + '.zip /awork05-1/'
sp.call(cpcp,shell=True)
