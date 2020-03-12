import pandas as pd
import os,sys,glob,psutil
import multiprocessing
#lists
rec_n = []
lig_n = []
ac = []
nac = []
sk = []
hh = []
rmsd = []
rfeats = []
rfeats1 = []
#dictionaries
rec_dic = {}
rec_dic2 = {}
lig_dic = {}
lig_dic2 = {}
cc_dic = {}
def reading_DB(x):
    hladb = pd.read_csv(x,sep='\t')
    indb = hladb[hladb['PDB ID'] == pdbid].dropna()
    hlaclass = indb['HLA CLASS']
    hlatype = indb['HLA TYPE']
    o_seq = indb['Peptide Sequence']
    seqlen = len(o_seq)
    hla = hlaclass.split('-')[1] + ''.join(hlatype.split(':'))

def native(x):
    recl = []
    pepl = []
    rec_ser = []
    lig_ser = []
    rec_n = []
    lig_n = []
    with open(RECLIB + hla + '/' + x + '_A.pdb','r') as R:
        for line in R.readlines():
            rec_n.append(line[7:13].strip())
            recl.append(line.strip())
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

    with open(PEPLIB + hla + '/' + x + '_B.pdb','r') as P:
        for line in P.readlines():
            lig_n.append(line[7:13].strip())
            pepl.append(line.strip())
            if line[13:17].strip() == 'OXT':
                pass
            else:
                if nnn == int(line[22:26].strip()):
                    lig_ser.append(line[13:17].strip())
                    lig_dic[nnn] = lig_ser
                else:
                    lig_dic[nnn] = lig_ser
                    nnn = int(line[22:26].strip())
                    lig_ser = []
                    lig_ser.append(line[13:17].strip())

    with open(x + '_native.pdb','w') as W:
        for line in recl:
            W.write(line + '\n')
        W.write('TER\n')
        for line in pepl:
            W.write(line + '\n')
        W.write('END\n')

def mutate(pdbid,inseq):
    mtd_lines = []
    os.environ['MMTSB'] = '/lwork01/mmtsb'
    mmtsb = os.environ['MMTSB']
    os.environ['PATH'] +=  os.environ['PATH'] + ':' + mmtsb + '/perl:' + mmtsb + '/bin'
    os.system('mutate.pl -seq 1:' + inseq + ' ' + pdbid + '_B.pdb > ' + inseq + '_mtd.pdb')
    with open(inseq + '_mtd.pdb','r') as MP:
        for line in MP.readlines():
            if not line.startswith('TER\n' or 'END\n') :
                mtd_lines.append(line[:21] + 'B ' + line[23:56] + '1' + line[57:])
            else:
                mtd_lines.append(line)
    del mtd_lines[-1]
    with open(pdbid + '_' + inseq + '_pep.pdb','w') as W:
        W.write('TER\n')
        for line in mtd_lines:
            W.write(line)
    os.system('cat ' + pdbid + '_A.pdb ' + pdbid + '_' + inseq + '_pep.pdb > ' + pdbid + '_' + inseq + '_inp.pdb')

def write_exec(pdb,seq):
    prepack_lines=['-s ' + pdb + '_' + seq + '_inp.pdb',\
                   '-out:path:all ' + RES,\
                   '-out:prefix PPK_',\
                   '-out:file:scorefile score_ppk_' + pdb + '_' + seq + '.sc',\
                   '-ex1',\
                   '-ex2aro',\
                   '-user_input_sc',\
                   '-flexpep_prepack',\
                   '-nstruct 1']
    firststep_lines =['-s ' + RES + '/PPK_' + pdb + '_' + seq + '_inp_0001.pdb',\
                  '-native ' + ROOT + 'Native/' + HLAclaty + '_native/' + pdb + '_native.pdb',\
                  '-out:path:all ' + RES + '/PDB_1ST',\
                  '-out:file:scorefile score_1st_' + pdb + '_' + seq + '.sc',\
                  'nstruct ' + nout,\
                  '-flexPepDocking:pep_refine',\
                  '-flexPepDocking:flexpep_score_only',\
                  '-ex1',\
                  '-ex2aro']
    secondstep_lines=['-s ' + RES + '/PPK_' + pdb + '_' + seq + '_inp_0001.pdb',\
                  '-native ' + ROOT + 'Native/' + HLAclaty + '_native/' + pdb + '_native.pdb',\
                  '-out:path:all ' + RES + '/PDB_2ND',\
                  '-out:file:scorefile score_2nd_' + pdb + '_' + seq + '.sc',\
                  'nstruct ' + nout,\
                  '-flexPepDocking:pep_refine',\
                  '-flexPepDocking:flexpep_score_only',\
                  '-ex1',\
                  '-ex2aro']
    with open('prepack_flags','w') as W:
        for line in prepack_lines:
            W.write(line + '\n')
    with open('run_flag1','w') as W:
        for line in firststep_lines:
            W.write(line + '\n')
    with open('run_flag2','w') as W:
        for line in secondstep_lines:
            W.write(line + '\n')
def prepare_input(x):
    if inseq != o_seq:
        mutate(x,inseq)
        write_exec(pdbid,inseq)
        RES = pdb + '_' + inseq + '_' +  hla
    else:
        os.system('cat ' + x + '_A.pdb ' +  x + '_B.pdb > ' + x + '_' + o_seq + '_inp.pdb')
        write_exec(x,o_seq)
        RES = pdb + '_' + o_seq + '_' + hla

def enva_working(x):
    inf = x.split('.')[0]
    os.system('csplit -f "%s_" %s.pdb \'/TER/\';sed \'s/ATOM  /HETATM/\' %s_01 > %s_02; cat %s_00 %s_02 > %s_het.pdb'%(inf,inf,inf,inf,inf,inf,inf))
    os.system('cp %s.pdb %s.pdb'%(inf,inf.split('_')[-3]))
    if not os.path.exists(inf + '_a.out'):
        os.system(GEAR + '/enva_rec.v1.1 -a ' + inf + '.pdb > ' + inf + '_a.out')
    if not os.path.exists(inf + '_b.out'):
        os.system(GEAR + '/enva_rec.v1.1 -b ' + inf + '_het.pdb > ' + inf + '_b.out')
    if not os.path.exists(inf + '_m.out'):
        os.system(GEAR + '/enva_rec.v1.1 -m ' + inf.split('_')[-3] + '.pdb ' + inf.split('_')[-3] + 'B > ' + inf + '_m.out')
    if not os.path.exists('*.env'):
        os.system(GEAR + '/enva_rec.v1.1 -e ' + inf + '.pdb B')

def atom_revise(x):
    rnn = 1
    nnn = 1
    olig_ser = {}
    orec_ser = {}
    lig_lines = []
    out_lines = []
    with open(x,'r') as F:
        for i in F.readlines():
            if i[21:23].strip() == 'A':
                if rnn == int(i[23:26].strip()):
                    orec_ser[i[13:17].strip()] = i
                    rec_dic2[rnn] = rec_ser
                else:
                    rec_dic2[rnn] = rec_ser
                    rnn = int(i[23:26].strip())
                    orec_ser = {}
                    orec_ser[i[13:17].strip()] = i
            elif i[21:23].strip() == 'B':
                if nnn == int(i[23:26].strip()):
                    olig_ser[i[13:17].strip()] = i
                    lig_dic2[nnn] = lig_ser
                else:
                    lig_dic2[nnn] = lig_ser
                    nnn = int(i[23:26].strip())
                    olig_ser = {}
                    olig_ser[i[13:17].strip()] = i
    feature = range(1,len(lig_dic.keys()) + 1)
    for i in feature:
        for j in lig_dic[i]:
            lig_lines.append(lig_dic2[i][j].strip())
    for i in rec_dic:
        for j in rec_dic[i]:
            out_lines.append(rec_dic2[i][j].strip())

    with open(x.split('.')[0] + '_rev.pdb','w') as W:
        for i,j in zip(out_lines,rec_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:] + '\n')
        W.write('TER\n')
        for i,j in zip(lig_lines,lig_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:] + '\n')
        W.write('END\n')



def txt_writing(x):
    ff_list = []
    hh_list = []
    tt_acc = []
    tt_phi = []
    tt_psi = []
    mat_lst = []
    ave_acc = []
    tot_lig_acc = 0
    ave_lig_acc = 0
    tot_lig_num = 0
    pdbl = []
    tttt = []
    newl = []

    for i in sorted(glob.glob('*rev.pdb')):
        os.system('csplit -f \'' + i + '_\' ' + i + ' \'/TER/\'')
        os.system(
            '/awork06-1/YKLee/c_script/rmsd_total ' + PEPLIB + '/' + inpdb + '_pep.pdb ' + i + '_01 bb >> ' + x + '_rmsd.txt')
    for f, g in zip(sorted(glob.glob('*a.out')), sorted(glob.glob('*m.out'))):
        sidx = []
        pdb = '_'.join([f.split('.')[0].split('_')[1], str(x),
                        f.split('.')[0].split('_')[5]])
        pdbl.append(pdb)
        adf = pd.read_csv(f, sep='\s+', skiprows=[0, 0], header=None)
        mdf = pd.read_csv(g, sep='\s+', header=None)
        mdf_filter = mdf[mdf[17] == 1].iloc[:, 9:16]
        for i in mdf_filter[9]:
            tot_lig_acc = tot_lig_acc + float(i)
            tot_lig_num += 1
            if tot_lig_num > 0:
                ave_lig_acc = tot_lig_acc / float(tot_lig_num)
            else:
                ave_lig_acc = 0
        ave_acc.append(ave_lig_acc)

        new_acc = {}
        num = 0
        for i in rfeats:
            new_acc['AA_' + cc_dic[str(i)]] = sum(
                [adf[adf[1] == int(i.split('_')[0])][adf[2] == i.split('_')[2]][adf[3] == i.split('_')[1]][11].values,
                 1])
        for i in mdf_filter.T:
            tidx = str(mdf_filter.T[i][13]) + '_' + str(mdf_filter.T[i][14]) + '_' + str(mdf_filter.T[i][15])
            if tidx in rfeats:
                tttt.append(tidx)
                zdx = rfeats.index(tidx) + 1
                sidx.append(zdx)
                num += 1
        newl.append(pd.DataFrame(new_acc))
        if not os.path.exists(pdb.split('.')[0] + '.ser'):
            with open(pdb.split('.')[0] + '.ser', 'w') as W:
                for sid in list(set(sidx)):
                    W.write(str(sid) + '\n')
        iskew_cmd = GEAR + '/iskew ' + pdb.split('.')[0] + '.ser >> total_' + x + '_sk.txt'
        subprocess.call(iskew_cmd, shell=True)
        ratio = float(num) / float(len(rfeats))
        mat_lst.append(ratio)
    new_ac = pd.concat(newl)
    new_ac['total_rec_acc'] = new_ac.sum(axis=1)
    new_ac['ave_lig_acc'] = ave_acc
    new_ac['%Match'] = mat_lst
    new_ac['PDB'] = pdbl
    new_ac.set_index('PDB').reset_index()

    new_ac.to_csv(x + '_nac.txt', sep='\t', index=False)
    nn_list = []
    for f in sorted(glob.glob('*b.out')):
        with open(f, 'r') as F:
            n = 0
            for line in F.readlines():
                tt = line[7:13].strip() + '_' + line[17:21].strip() + '_' + line[13:17].strip()
                if tt in rfeats:
                    n += 1
            nn_list.append(n)
        with open(f, 'r') as F:
            hh_list.append(len(F.readlines()))
        ff_list.append(f.split('.')[0].split('_')[1] + '_' + str(x) + '_' + f.split('.')[0].split('_')[-4])
    with open(x + '_hh.txt', 'w') as W:
        W.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
        for i, j, k in zip(ff_list, hh_list, nn_list):
            W.write(str(i) + '\t' + str(j) + '\t' + str(k) + '\n')
    acc_columns = ['P%d' % i for i in range(1, seqlen + 1)]
    phi_columns = ['PHI%d' % i for i in range(1, seqlen + 1)]
    psi_columns = ['PSI%d' % i for i in range(1, seqlen + 1)]
    pdbl = []
    for f in sorted(glob.glob('*.env')):
        pdb = '_'.join([sam.split('_')[1], str(x), f.split('.')[0].split('_')[5]])
        pdbl.append(pdb)
        with open(f, 'r') as F:
            acc = []
            psi = []
            phi = []
            for line in F.readlines():
                if line.find('chain') and line.startswith('ATOM'):
                    envs = line.split()[9:]
                    acc.append(envs[2])
                    phi.append(envs[4])
                    psi.append(envs[5])
            tt_acc.append(pd.DataFrame(acc, index=acc_columns).T)
            tt_phi.append(pd.DataFrame(phi, index=phi_columns).T)
            tt_psi.append(pd.DataFrame(psi, index=psi_columns).T)
    pdbdf = pd.DataFrame(pdbl, columns=['PDB'])
    acc_df = pd.concat(tt_acc, ignore_index=True).reindex()
    phi_df = pd.concat(tt_phi, ignore_index=True).reindex()
    psi_df = pd.concat(tt_psi, ignore_index=True).reindex()
    act_df = pd.concat([pdbdf, acc_df, phi_df, psi_df], axis=1)
    act_df.to_csv(x + '_ac_ct.txt', sep='\t', index=False)
def reT(x):
    Olines = []
    recl = []
    pepl = []
    with open(x,'r') as F:
        for line in F.readlines():
            tline = line.strip()
            if tline.startswith('ATOM'):
                if tline[-1] != 'H':
                    Olines.append(tline)
    for line in Olines:
        if line[21] == 'A':
            recl.append(line)
        elif line[21] != 'A':
            pepl.append(line[:21] + 'B ' + line[23:])
    with open(x.split('.')[0] + 'red.pdb','w') as W:
        for line in recl:
            W.write(line + '\n')
        W.write('TER\n')
        for line in pepl:
            W.write(line + '\n')
def enva_part(y):
    os.chdir(y)
    os.system('rm -rf *rev* *red* *.ser *.txt')
    for i in sorted(glob.glob('*[0-9].pdb')):
        reT(i)
        atom_revise(i)
    list1 = sorted(glob.glob('*rev.pdb'))
    pool1 = mutlprocessing.Pool(ncpu)
    pool1.map(enva_working(list1))
    pool1.close()
    pool1.join()
    txt_writing(y.split('_')[1])
    ac.append(pd.read_csv(y.split('_')[1] + '_ac_ct.txt', sep='\t'))
    nac.append(pd.read_csv(y.split('_')[1] + '_nac.txt', sep='\t'))
    sk.append(pd.read_csv('total_' + y.split('_')[1] + '_sk.txt', sep='\t', header=None).fillna('1.000000'))
    hh.append(pd.read_csv(y.split('_')[1] + '_hh.txt', sep='\t'))
    rmsd.append(pd.read_csv(y.split('_')[1] + '_rmsd.txt', sep='\t', header=None))
    os.chdir('../')



pdbid = sys.argv[1]
inseq = sys.argv[2]
nout = sys.argv[3]
dbfile = sys.argv[4]
reading_DB(dbfile)
ncpu = str(psutil.cpu_count()-1)
feats = pd.read_csv(PDBLIB + '/' + pdbid + '.out',sep='\t',header=None)
feats_filter = feats[feats[17] == 1].iloc[:,13:16]
for i in feats_filter.values.tolist():
    rfeats.append(str(i[0]) + '_' + '_'.join(i[1:]))
feats1_filter = feats[feats[17] ==1].iloc[:,12:16]
for i in feats1_filter.values.tolist():
    rfeats1.append(str(i[0]) + '_' + '_'.join(i[2:]))
for i,j in zip(rfeats,rfeats1):
    cc_dic[i] = j

if int(ncpu) >= 17 :
    os.environ['rosetta'] = '/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin/'
    os.environ['rosetta_clean'] = '/lwork01/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/'
    os.environ['neogear'] = '/lwork01/neoscan_gear/'
    gear = os.environ['neogear']
    rosetta = os.environ['rosetta']
    rosetta_clean = os.environ['rosetta_clean']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta + ':' + os.environ['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta_clean + ':' + os.envirion['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + gear + ':' + os.environ['PATH']
    os.system('scp -rpq user1@10.1.5.6:/awork06-1/YKLee/pdpdb/Neoscan_V2/Native /lwork01/')
    os.system('scp -rpq user1@10.1.5.6:/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/receptor/ /lwork01/')
    os.system('scp -rpq user1@10.1.5.6:/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/peptide/ /lwork01/')
    PDBLIB = '/lwork01/Native/'
    PEPLIB = '/lwork01/peptide/'
    RECLIB = '/lwork01/receptor/'
else:
    os.environ['rosetta'] = '/awork06-1/rosetta_src_2019.40.60963_bundle/main/source/bin/'
    os.environ['rosetta_clean'] = '/awork06-1/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/'
    os.environ['neogear'] = '/awork06-1/neoscan_gear/'
    gear = os.environ['neogear']
    rosetta = os.environ['rosetta']
    rosetta_clean = os.environ['rosetta_clean']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta + ':' + os.environ['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + rosetta_clean + ':' + os.envirion['PATH']
    os.environ['PATH'] += ':' + os.environ['PATH'] + gear + ':' + os.environ['PATH']
    PDBLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/Native'
    RECLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/receptor/'
    PEPLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/peptide/'
#rosetta clean py, divide by chain
os.system('cat ' + RECLIB + pdbid + '_rec.pdb ' + PEPLIB + pdbid + '_pep.pdb > ' + pdbid + '.pdb')
os.system('clean_pdb.py ' + pdbid + '.pdb A')
os.system('clean_pdb.py ' + pdbid + '.pdb B')
native(pdbid)
prepare_input(pdbid)

os.system('FlexPepDocking.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/preapack_flags > ' + RES + '/prepack.log')
for i,j in zip(['run_flag1','run_flag2'],['run1.log','run2.log']):
    os.system('mpirun -np ' + ncpu + ' FlexPepDocking.mpi.linuxgccrelease -database ' + ROSETTA_DB + ' @' + RES + '/' + i + '> ' + RES + j)

for nn in ['PDB_1ST','PDB_2ND']:
    enva_part(nn)
os.chdir(RES + '_energy_matrix/')
df_ac = pd.concat(ac)
df_nac = pd.concat(nac)
df_sk = pd.concat(sk)
df_hh = pd.concat(hh)
df_rmsd = pd.concat(rmsd)
df_ac.to_csv('total_ac.txt',sep='\t',index=False)
df_rmsd.to_csv('total_rmsd.txt',sep='\t',index=False)
df_nac.to_csv('total_nac.txt',sep='\t',index=False)
df_sk.to_csv('total_sk.txt',sep='\t',index=False)
df_hh.to_csv('total_hh.txt',sep='\t',index=False)
df_rmsd.columns = ['pdb','reference','rmsd']
total_df = [df_hh,df_nac,df_sk,df_ac]
df_final1 = reduce(lambda left,right: pd.merge(left,right,on=['PDB'],how='outer'),total_df)
df_final1['rmsd'] = df_rmsd['rmsd']
df_final2 = df_final1.set_index('rmsd').reset_index().set_index('PDB').reset_index()
df_final2.to_csv(RES + '_total.txt',sep='\t',index=False,na_rep='-')
os.chdir('../')