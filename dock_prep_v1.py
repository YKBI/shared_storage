import os
import sys
import subprocess
import shutil

seq = sys.argv[1]
HLAclass = sys.argv[2]
HLAtype = sys.argv[3]
HLAclaty = HLAclass[-1]+HLAtype
ROOT = '/awork06-1/YKLee/NAVS/'
REC_cwd = ROOT + 'receptor/' + HLAclaty + '/'
PEP_cwd = ROOT + 'peptide/' + HLAclaty + '/'
RES = '%s_%s_%s'%(seq,HLAclass,HLAtype)
try:
    if not os.path.exists('%s_%s_%s'%(seq,HLAclass,HLAtype)):
        os.mkdir('%s_%s_%s'%(seq,HLAclass,HLAtype),0777)
except OSError :
    pass

pdbs = []
with open(ROOT + HLAclass + '_list/' + HLAtype + '.list','r') as f:
    for line in f.readlines():
        pdbs.append(line.strip())

for pdb in pdbs:
    mutation = '/awork06-1/YKLee/pl_script/toolset/perl/mutate.pl -seq 1:' + seq + " " + PEP_cwd + pdb + '_pep.pdb > ' +  RES + '/' + pdb + '_' + seq + '_temp.pdb'
    subprocess.call(mutation,shell=True)
    p_file = RES + '/' + pdb + '_' + seq + '.pdb'
    with open(RES + '/' + pdb + '_' + seq + '_temp.pdb','r') as f1:
        OLines = ['TER\n']
        for line in f1.readlines():
            if not line.startswith('TER\n' or 'END\n'):
                oLine = line[:21] + 'B ' + line[23:56] + '1' + line[57:]
                OLines.append(oLine)
            else:
		oLine = line
		OLines.append(oLine)
    	del OLines[-1]
	OLines.append('END')
    with open(p_file,'w') as f2:
	for line in OLines:
		f2.write(line)
    cat = 'cat ' + REC_cwd +pdb + '_rec.pdb ' + p_file + ' > ' + RES + '/' + pdb + '_' + seq + '_inp.pdb'
    subprocess.call(cat,shell=True)

