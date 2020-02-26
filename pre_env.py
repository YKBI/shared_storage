
import os
import sys
import subprocess as sp

inp_dir = sys.argv[1]
HLA = '_'.join(inp_dir.split('_')[1:])
seq = inp_dir.split('_')[0]
HLAtype = ''.join(inp_dir.split('_')[1:]).split('-')[1]#''.join(inp_dir.split('_')[1:])
o_dir = '/home/user1/' + seq
print o_dir
pdb_dic = {'HLA-A_0101':['3bo8','4nqx','4nqv','1w72','5brz'],\
	'HLA-A_0201':['1oga','2gtw'],\
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

def tt_v3(seq,pdb,HLAtype,o_dir):
	cmd1 = 'time python /awork06-1/YKLee/py_script/tt_v3.py %s %s %s %s'%(seq,pdb,HLAtype,o_dir)
	print cmd1
	sp.call(cmd1,shell=True)

try:
	if not os.path.exists(o_dir):
		os.makedirs(o_dir)
except OSError:
	pass

os.chdir(inp_dir)
print os.getcwd()
for pdb in pdb_dic[HLA]:
	tt_v3(seq,pdb,HLAtype,o_dir)


