#/usr/bin/env python

import os
import sys
import shutil
import subprocess

def help():
	print "print help usage\n"
	print "Usage: python %s -s [ sample name ] -i [ input folder] -r [ ref version ] -f [ format ] -d [ cancer or non-cancer]\n"%sys.argv[0]
	return

if len(sys.argv) == 1:
	print "Supply suitable option\n"
	help()
else:
	sam = sys.argv[1]
	sam1 = sys.argv[2]

	try:
		if not os.path.exists(sam):
			os.mkdir(sam,0777)
	except OSError:
		pass

	wdir = os.getcwd()
	shutil.copy(sam + '.pdb',sam)
	os.chdir(sam)

	with open('rec.pdb','w') as recf:
		with open(sam + '.pdb','r') as comf:
			lines = comf.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
					recf.write(line)
		recf.write('END\n')

	with open('lig.pdb','w') as ligf:
		with open(sam + '.pdb','r') as comf:
			lines = comf.readlines()
			for line in lines:
				if line.startswith('HETATM') > 0 :
					ligf.write(line)
		ligf.write('END\n')

	cmd = 'prep_01_simulation_protein_small_molecule_complex_v1.4_min.py -ir rec.pdb -il lig.pdb -t %s -fix-gpu-id %s -auto'%(sam,sam1)
	subprocess.call(cmd,shell=True)

	os.chdir(wdir)
