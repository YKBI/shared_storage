import sys
import os
import subprocess as sp
f = sys.argv[1]
nuli = []
qwli = []
#fna = '_'.join(f.split('.'[0]))+"_B.pdb"
nudic = {}
Lines = []
with open(f,'r') as a:
	reads = a.readlines()
	for i in reads:
		if i.startswith('TER'):
			T_num = int(i[22:26].strip()) 
with open(f,'r') as b:
	for line in b.readlines():
		if line.startswith('ATOM'):		
			if T_num < int(line[22:26].strip()): nuli.append(int(line[22:26].strip()))

for i,j in zip(range(1,10),sorted(list(set(nuli)))):
	nudic[j] = i
print nudic
with open(f,'r') as c:
	for line in c.readlines():
		if line.startswith('ATOM'):
			if T_num < int(line[22:26].strip()): Lines.append(line[:23] + str(nudic[int(line[22:26].strip())]).rjust(3) + line[26:])
			else : Lines.append(line)
		elif line.startswith('TER'): Lines.append('TER\n')
with open(f+'new','w') as W:
	for eline in Lines:
		W.write(eline)
	W.write('END')
'''
rm = "rm -rf %s"%f
mv = "mv %s %s"%(f+'new',f)
sp.call(rm,shell=True)
sp.all(mv,shell=True)'''
