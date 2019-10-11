import sys
import os
import shutil as sh
import glob as gb


fs = gb.glob('*.pdb.*')
num1 = int(sys.argv[1])#helix part 1 start
num2 = int(sys.argv[2])#helix part 1 end
num3 = int(sys.argv[3])#helix part 2 start
num4 = int(sys.argv[4])#helix part 2 end
try:
	os.makedirs('%s_v8/traj_1/pdb_from_prod/'%'/'.join(os.getcwd().split('/')[:5]))

except:
	for f in fs:
		lines = []
		with open(f,'r') as F:
			reads = F.readlines()
			for i in reads:
				if i.startswith('ATOM' or 'TER' or 'END'):
					if num1 <= int(i[22:26].strip())<= num2:
						lines.append(i)
					elif num3 <= int(i[22:26].strip()) <= num4:
						lines.append(i)
					#elif num4+3< int(i[22:26].strip()):
					#	lines.append(i)
		for i in lines:
			with open('%s_new'%f,'a') as N:
				N.write(i)
filelist = gb.glob('*new')
for i in filelist:
	sh.move(i,'%s_v8/traj_1/pdb_from_prod/%s'%('/'.join(os.getcwd().split('/')[:5]),i.split('_new')[0]))
