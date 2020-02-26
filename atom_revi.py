
import sys
import os
pdbna = sys.argv[1]

rec_n = []
lig_n = []
fff = []
out_lines = []
het_lines = []
ddd = {}
het_dic = {}

#orig pdb reading
def reading(x):
    with open(x + '.pdb','r') as F:
        for i in F.readlines():
            if i.startswith('ATOM'):
                rec_n.append(int(i[22:26].strip()))
            elif i.startswith('HETATM'):
                lig_n.append(int(i[22:26].strip()))
                fff.append(i[13:17].strip())

#pdb of minimize reading
def reading2(x):
    with open(x + '_min_conv.pdb','r') as F:
        for i in F.readlines():
            if i.startswith('ATOM'):
                out_lines.append(i[:22] + str(ddd[int(i[22:26].strip())]).rjust(4) + i[26:])
            if i.startswith('HETATM'):
                het_lines.append(i[:22] + str(ddd[int(i[22:26].strip())]).rjust(4) + i[26:])
        out_lines.append('TER\n')

#revise pdb writing
def writing(x):
    with open('../' + x + '_rev.pdb','w') as W:
        for i in out_lines:
            W.write(i)
        for j in het_lines:
            W.write(j[:13] + het_dic[j[13:17].strip()].ljust(4) + j[17:])
        W.write('END\n')

#activation
os.chdir(pdbna)

reading(pdbna)
rec_nn = sorted(list(set(rec_n))) + list(set(lig_n))
for i,j in zip(range(1,len(rec_nn)+1),rec_nn):
    ddd[i] = j
reading2(pdbna)
for i,j in zip(het_lines,fff):
    het_dic[i[13:17].strip()] = j
writing(pdbna)
os.chdir('../')
