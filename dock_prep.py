#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import string


sam = sys.argv[1] # seq
sam1 = sys.argv[2] # HLA class
sam2 = sys.argv[3] # HLA type
RES = sam + "_" + sam1 + "_" + sam2
ROOT = '/awork04-2/NAVS/'

try:
	if not os.path.exists(sam + "_" + sam1 + "_" + sam2):
		os.mkdir(sam + "_" + sam1 + "_" + sam2,0777)
except OSError :
	pass

pdbs =[]
with open(ROOT + sam1 + "_list/" + sam2 + ".list",'r') as f:
	lines = f.readlines()
	for line in lines:
		pdbs.append(line[:-1])

for pdb in pdbs:
	cmd = "mutate.pl -seq 1:" + sam + " " + ROOT + "/" + sam1 + "/" + sam2 + "/" + pdb + "/" + pdb + "_pep.pdb > " + RES + "/" + pdb + "_" + sam + "_temp.pdb"
#	print cmd
	subprocess.call(cmd,shell=True)
	pep_file = RES + "/" + pdb + "_" + sam + ".pdb"
	with open(pep_file,'a') as f2:
		f2.write("TER\n")
		with open(RES + "/" + pdb + "_" + sam + "_temp.pdb",'r') as f1:
			lines = f1.readlines()
			for line in lines :
				if line.find('TER')<0 :
					if line.find('END')<0 :
						TEXT = line[:21] + 'B ' + line[23:] 
						f2.write(TEXT)
					else:
						f2.write(line[:-1])
	cmd1 = "cat " + ROOT + "/" + sam1 + "/" + sam2 + "/" + pdb + "/" + pdb + "_rec.pdb " + pep_file + " > " + RES + "/" + pdb + "_" + sam + "_inp.pdb"
	subprocess.call(cmd1,shell=True)
