#/usr/bin/env python

import os
import sys
import getopt
import subprocess
import shutil
import string
import pandas as pd
from pandas import DataFrame
from collections import Counter
import matplotlib.pyplot as plt
import math
#from sklearn import datasets, linear_model
#import matplotlib.pyplot as plt
#import statsmodels.api as sm
import numpy as np
import glob

sam = sys.argv[1]
sam1 = sys.argv[2]

df = pd.read_csv(sam,sep='\t')
serial = df.index.tolist()

lig = sam.split('.')[0].split('_')[1]
comp = sam.split('.')[0]

feats = []
with open(sam1,'r') as ff :
	lines = ff.readlines()
	for line in lines:
		feats.append(line[:-1])

#df1 = df[df['rmsd']<=1.0]

for feat in feats:
	plt.clf()
	plt.scatter(serial,df[feat],label=lig)
	plt.xlabel('serial')
	plt.ylabel(feat)
	plt.legend()
	plt.savefig(comp + '_' + feat + '.png')
