import sys
import os

f = sys.argv[1]
num1 = int(sys.argv[2])
num2 = int(sys.argv[3])
num3 = int(sys.argv[4])
num4 = int(sys.argv[5])
helix1 = []
helix2 = []
with open(f,'r') as F:
	reads = F.readlines()
	for i in reads:
		if num1 <= int(i[23:27].strip()) <= num2:
			helix1.append(i[23:27].strip())
		elif num3 <= int(i[23:27].strip()) <= num4:
			helix2.append(i[23:27].strip())
helix1_v1 = list(set(helix1))
helix2_v1 = list(set(helix2))
helix1_v1.sort()
helix2_v1.sort()
features = helix1_v1 + helix2_v1
os.chdir('/awork06-1/YKLee/NAVS/feat_list/')
with open('%s.list'%f.split('/')[-1].split('_rec')[0],'a') as A:
	for i in features:
		A.write(i+'\n')
