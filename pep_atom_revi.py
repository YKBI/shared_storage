
import sys
import os
#pdbna = sys.argv[1]
samin1 = sys.argv[1]
samin2 = sys.argv[2]
rec_dic = {}
rec_n = []
fff = []
out_lines = []
het_lines = []
ddd = {}
het_dic = {}
het_dic2 = {}
rec_dic2 = {}
rec_n = []
lig_n = []
feature = []
feature_dic = {}
#orig pdb reading
def reading(x):
    rec_ser = []
    lig_ser = []
    nnn = 1
    rnn = 1
    with open('../' + x + '_A.pdb','r') as F: #'/awork06-1/YKLee/pdpdb/Neoscan_V2/receptor_pdb/' + x + '_A.pdb','r') as F:
        for i in F.readlines():
	    if i[21:23].strip() == 'A':
		if i[13:17].strip() == 'OXT':
		    pass
		else:
                    if rnn == int(i[22:26].strip()):
		        rec_ser.append(i[13:17].strip())
		        rec_dic[rnn] = rec_ser
		    else:
		        rnn = int(i[22:26].strip())
		        rec_ser = []
		        rec_ser.append(i[13:17].strip())
    with open('../' + x + '_B.pdb','r') as F : #'/awork06-1/YKLee/pdpdb/Neoscan_V2/peptide_pdb_conv/' + x + '_B.pdb','r') as F:
	for i in F.readlines():
	    if i[21:23].strip() == 'B':
		if nnn == int(i[22:26].strip()):
		    lig_ser.append(i[13:17].strip())
		    het_dic[nnn] = lig_ser
		else:
		    #print nnn
		    het_dic[nnn] = lig_ser
                    nnn = int(i[22:26].strip())
		    lig_ser = []
		    lig_ser.append(i[13:17].strip())
	        
def origin(x):
    rec_ser = []
    #rec_n = []
    lig_ser = []
    #lig_n = []
    nnn =1
    rnn =1
    with open(x,'r') as F:
	for line in F.readlines():
	    if line[21:23].strip() == 'A':
		rec_n.append(line[7:13].strip())
		if line[13:17].strip() == 'OXT':
		    pass
		else:
		    if rnn == int(line[22:26].strip()):
			rec_ser.append(line[13:17].strip())
			rec_dic[rnn] = rec_ser
		    else:
			rnn = int(line[22:26].strip())
			rec_ser = []
			rec_ser.append(line[13:17].strip())
	    elif line[21:23].strip() == 'B':
		lig_n.append(line[7:13].strip())
		if line[13:17].strip() == 'OXT':
		    pass
		else:
		    if nnn == int(line[22:26].strip()):
			lig_ser.append(line[13:17].strip())
			het_dic[nnn] = lig_ser
		    else:
			het_dic[nnn] = lig_ser
			nnn = int(line[22:26].strip())
			lig_ser = []
			lig_ser.append(line[13:17].strip())
#pdb of minimize reading
def reading2(x):
    
    rnn = 1
    nnn = 1
    lig_ser = {}
    rec_ser = {}
    #with open(x + '_min_red.pdbnew','r') as F:
    with open(x,'r') as F:
        for i in F.readlines():
            if i[21:23].strip() =='A': #.startswith('ATOM'):
        	if rnn == int(i[23:26].strip()):
		    rec_ser[i[13:17].strip()] = i
		    rec_dic2[rnn] = rec_ser
		else:
		    rec_dic2[rnn] = rec_ser
		    rnn = int(i[23:26].strip())
		    rec_ser = {}
		    rec_ser[i[13:17].strip()] = i
            elif i[21:23].strip() == 'B': #.startswith('HETATM'):
	        if nnn == int(i[23:26].strip()):
                    lig_ser[i[13:17].strip()] = i
		    het_dic2[nnn] = lig_ser
		else:
		    het_dic2[nnn] = lig_ser
		    nnn = int(i[23:26].strip())
		    lig_ser = {}
		    lig_ser[i[13:17].strip()] = i
def reading3(x):
    rnn = 1
    nnn = 1
    lig_ser = {}
    rec_ser = {}
    with open(x + '_00','r') as F:
	for i in F.readlines():
	    if i.startswith('ATOM'):
		if i[13:17].strip() == 'OXT':
		    pass
		else:
	            if rnn == int(i[23:26].strip()):
		        rec_ser[i[13:17].strip()] = i
		        rec_dic2[rnn] = rec_ser
	            else:
		        rec_dic2[rnn] = rec_ser
		        rnn = int(i[23:26].strip())
	 	        rec_ser = {}
		        rec_ser[i[13:17].strip()] = i
    with open(x + '_01','r') as F:
	for i in F.readlines():
	    if i.startswith('ATOM'):
		if i[13:17].strip() == 'OXT':
		    pass
		else:
	            if nnn == int(i[23:26].strip()):
		        lig_ser[i[13:17].strip()] = i
		        het_dic2[nnn] = lig_ser
	            else:
		        het_dic2[nnn] = lig_ser
		        nnn = int(i[23:26].strip())
		        lig_ser = {}
		        lig_ser[i[13:17].strip()] = i
#revise pdb writing
def writing1(x):
    with open('_'.join(x.split('_')[:2]) + '_rev.pdb','w') as W:
        for i,j in zip(out_lines,rec_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:21] + 'A ' + i[23:] + '\n')
	W.write('TER\n')
        for i,j in zip(het_lines,lig_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:21] + 'B ' + i[23:] + '\n') #j[:13] + het_dic[j[13:17].strip()].ljust(4) + j[17:])
        W.write('END\n')
def writing2(x):
    with open(x.split('.')[0] + '_rev.pdb','w') as W:
	for i,j in zip(out_lines,rec_n):
	    W.write(i[:7] + j.rjust(4) + '  ' + i[13:] + '\n')
	W.write('TER\n')
	for i,j in zip(het_lines,lig_n):
	    W.write(i[:7] + j.rjust(4) + '  ' + i[13:] + '\n')
	W.write('END\n') 
#os.system('csplit -f \'%s_\' %s \"/TER/\"'%(samin2,samin2))
#activation


if os.path.exists('feat.list'):
    with open('feat.list','r') as F:
	for line in F.readlines():
	    feature.append(line.strip())
    origin(samin1)
    reading3(samin2)
    for i,j in zip(het_dic.keys(),feature):
	feature_dic[j] = i
    for i in feature:
	for j in het_dic[feature_dic[i]]:
	    het_lines.append(het_dic2[int(i)][j].strip())

else:
    #feature = range(1,len(het_dic.keys())+1)
    origin(samin1)
    reading2(samin2)
    feature = range(1,len(het_dic.keys())+1)
    for i in feature:
	for j in het_dic[i]:
	    het_lines.append(het_dic2[i][j].strip())
'''
origin(samin1)
reading3(samin2)
feature = range(1,len(het_dic.keys())+1)

for i in feature:
	for j in het_dic[i]:
		het_lines.append(het_dic2[i][j].strip())
'''
for i in rec_dic:
    for j in rec_dic[i]:
	out_lines.append(rec_dic2[i][j].strip())

if os.path.exists('feat.list'):
    writing1(samin2)
else:
    writing2(samin2)

writing2(samin2)
#os.chdir('../')'''
