#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import string
import glob

def remove_duplicates(li):
	my_set = set()
	res = []
	for e in li:
		if e not in my_set:
			res.append(e)
			my_set.add(e)
	return res

def help():
	print "Usage: python %s [ seq ] [ PDB ] [ HLA-type ] [ final output folder]\n"%sys.argv[0]
	return

if len(sys.argv) == 1:
	help()
else :
	GEAR = '/awork06-1/neoscan_gear'
	sam = sys.argv[1]  # Sequence
	sam1 = sys.argv[2] # PDB
	sam2 = sys.argv[3] # HLA-type
	sam3 = sys.argv[4] # final output folder

	try:
		if not os.path.exists(sam3):
			os.mkdir(sam3,0777)
	except OSError:
		pass

	try:
		if not os.path.exists(sam3 + '/orig_pdb'):
			os.mkdir(sam3+ '/orig_pdb',0777)
	except OSError:
		pass

	try:
		if not os.path.exists(sam3 + '/score'):
			os.mkdir(sam3+ '/score',0777)
	except OSError:
		pass

#	try:
#		if not os.path.exists(sam4 + '/' + sam3):
#			os.mkdir(sam4+ '/' + sam3,0777)
#	except OSEkkrror:
#		pass

	wdir = os.getcwd()
	os.chdir( sam + '_' + sam1 + '_' + sam2 + '/PDB_1ST')
	fpdbs = []
	fpdbs = glob.glob('*[0-9].pdb')
	for fpdb in fpdbs :
		ser = int(fpdb.split('.')[0].split('_')[5])
		nfpdb = sam1 + '_1ST_' + str(ser)
		try :
			if not os.path.exists(sam3 + '/' + nfpdb):
				os.mkdir(sam3 + '/' + nfpdb,0777)
		except OSError:
			pass
		idx = 0 
		with open(sam3 + '/orig_pdb/' + nfpdb + '.pdb','a') as f :
			with open(fpdb,'r') as f1 :
				lines = f1.readlines()
				for line in lines :
					if line.startswith('ATOM') > 0  and line.find('OXT') < 0 :
						if idx == 0 and line[21]!='A' :
							TEXT = line[:21] + 'A' + line[22:]
							f.write(TEXT)
						elif idx > 0 and line[21]!='B' :
							TEXT = line[:21] + 'B' + line[22:]
							f.write(TEXT)
						else :
							f.write(line)
					elif line.startswith('TER') > 0 :
						if idx == 0 :
							f.write(line)
						idx = idx + 1 
			f.write('END\n')

		cmd = GEAR + '/reduce -Trim ' + sam3 + '/orig_pdb/' + nfpdb + '.pdb > ' + sam3 + '/orig_pdb/' + nfpdb + '_red.pdb'
		subprocess.call(cmd,shell=True)

		with open(sam3 + '/' + nfpdb +'/' + nfpdb + '_rec1.pdb','a') as rf :
			with open(sam3 + '/orig_pdb/' + nfpdb + '_red.pdb','r') as of :
				lines = of.readlines()
				for line in lines :
					if line.startswith('ATOM') > 0 and line[21]=='A' :
						rf.write(line)
			rf.write('TER\n')

		with open(sam3 + '/' + nfpdb + '/' + nfpdb + '_pep.pdb','a') as pf :
			with open(sam3 + '/orig_pdb/' + nfpdb + '_red.pdb','r') as of :
				lines = of.readlines()
				for line in lines :
					if line.startswith('ATOM') > 0 and line[21]=='B':
						pf.write(line)								
			pf.write('END\n')

	os.chdir('../PDB_2ND')
	lpdbs = []
	lpdbs = glob.glob('*[0-9].pdb')
	for lpdb in lpdbs :
		ser = int(lpdb.split('.')[0].split('_')[5])
		nfpdb = sam1 + '_2ND_' + str(ser)
		try :
			if not os.path.exists(sam3 + '/' + nfpdb):
				os.mkdir(sam3 + '/' + nfpdb,0777)
		except OSError:
			pass
		idx = 0 
		with open(sam3 + '/orig_pdb/' + nfpdb + '.pdb','a') as f :
			with open(lpdb,'r') as f1 :	
				lines = f1.readlines()
				for line in lines :
					if line.startswith('ATOM') > 0  and line.find('OXT') < 0 :
						if idx == 0 and line[21]!='A' :
							TEXT = line[:21] + 'A' + line[22:]
							f.write(TEXT)
						elif idx > 0 and line[21]!='B':
							TEXT = line[:21] + 'B' + line[22:]	
							f.write(TEXT)
						else :
							f.write(line)	
					elif line.startswith('TER') > 0 :
						if idx == 0 :
							f.write(line)
						idx = idx + 1
			f.write('END\n')
	
		cmd = GEAR + '/reduce -Trim ' + sam3 + '/orig_pdb/' + nfpdb + '.pdb > ' + sam3 + '/orig_pdb/' + nfpdb + '_red.pdb'
		subprocess.call(cmd,shell=True)

		with open(sam3 + '/' + nfpdb + '/' + nfpdb + '_rec1.pdb','a') as rf :
			with open(sam3 + '/orig_pdb/' + nfpdb + '_red.pdb','r') as of :
				lines = of.readlines()
				for line in lines :
					if line.startswith('ATOM') > 0 and line[21]=='A' :
						rf.write(line)
			rf.write('TER\n')

		with open(sam3 + '/' + nfpdb + '/' + nfpdb + '_pep.pdb','a') as pf :	
			with open(sam3 + '/orig_pdb/' + nfpdb + '_red.pdb','r') as of :
				lines = of.readlines()
				for line in lines:
					if line.startswith('ATOM') > 0 and line[21]=='B':
						pf.write(line)		
			pf.write('END\n')	


	os.chdir('../')
	with open(sam3 + '/score/' + sam1 + '_total_score.tsv','a') as sf :
		with open('total_score.tsv','r') as tf :
			lines = tf.readlines()
			for line in lines :
				if line.startswith('description') > 0 :
					sf.write('PDB\trmsBB\ttotal_score\n')
				else :
					cols = line.split('\t')
					order = cols[0].split('/')[0].split('_')[1]
					ser = int(cols[0].split('/')[1].split('.')[0].split('_')[5])
					fname = sam1 + '_' + order + '_' + str(ser)
					sf.write('%s\t%s\t%s'%(fname,cols[1],cols[2]))
	os.chdir(wdir)
