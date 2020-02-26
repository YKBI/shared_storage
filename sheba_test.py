import os
import sys
import pandas
import multiprocessing
import glob
from functools import partial
import pandas as pd

def sheba(a,b):
	#outform1 = '.'.join(['.'.join(a.split('.')[:2]),a.split('.')[3],a.split('.')[2]])
	outform2 = '.'.join(a.split('.')[:3]) 
	#os.system('cp ' + a + ' ' + outform1)
	os.system(GEAR + '/sheba_01 -x ' + b + ' ' + a) #PDBLIB + '/' + inpdb + '_native.pdb ' + outform1)
	os.system(GEAR + '/sheba_01 -t ' + '.'.join(a.split('.')[:3])+ '.trf ' + a)
	os.system('mv ' + a + '.pdb ' + outform2 + '_tr.pdb')
	os.system('csplit -f \'%s_\' %s.pdb \'/TER/\''%(outform2,outform2))

def enva_working(pdb):
        trpdb = '.'.join(pdb.split('.')[:3])#'_'.join(pdb.split('_')[:2])
	#os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py %s_rev.pdb A'%trpdb)
	os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py %s.pdb B'%trpdb)
	os.system('csplit -f \'%s_\' %s.pdb \'/TER/\''%(trpdb,trpdb))
	os.system('cat ../../%s_rec.pdb %s_B.pdb > %s_new.pdb'%(trpdb.split('.')[0],trpdb,trpdb))#(HLAclass,trpdb.split('.')[0],trpdb,trpdb))
        os.system('sed \'s/ATOM  /HETATM/\' %s_01 > %s_02 ;cat /awork06-1/YKLee/pdpdb/Neoscan_V2/v2/revise_pdb/receptor/%s/%s_rec.pdb %s_02 > %s_het.pdb'%(trpdb,trpdb,HLAclass,trpdb.split('.')[0],trpdb,trpdb))
        
	#os.system('cp %s.pdb %s.pdb'%(trpdb,trpdb.split('.')[2]))
        if not os.path.exists(inpdb + '_a.out'):
                os.system(GEAR + '/enva_rec.v1.1 -a ' + trpdb + '_new.pdb > ' + trpdb + '_a.out')
        if not os.path.exists(inpdb + '_b.out'):
                os.system(GEAR + '/enva_rec.v1.1 -b ' + trpdb + '_het.pdb > ' + trpdb + '_b.out')
        if not os.path.exists(inpdb + '_m.out'):
                os.system(GEAR + '/enva_rec.v1.1 -m ' + trpdb + '_new.pdb ' + trpdb.split('.')[0] + 'B > ' + trpdb + '_m.out')
	if not os.path.exists('*.env'):
		os.system(GEAR + '/enva_rec.v1.1 -e ' + trpdb + '_new.pdb B')
        #os.system('rm -rf ' + inpdb.split('_')[-3] + '.pdb')
	
def original_enva(pdb):
	trpdb = '.'.join(pdb.split('.')[:3]) #'_'.join(pdb.split('_')[:2])
	os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py %s.pdb A'%trpdb)
        os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py %s.pdb B'%trpdb)
	
	#os.system('cat %s_db %s_01 > %s_new.pdb'%(trpdb,trpdb,trpdb))
	os.system('sed \'s/ATOM  /HETATM/\' %s_B.pdb > %s_02 ;cat %s_A.pdb %s_02 > %s_het.pdb'%(trpdb,trpdb,trpdb,trpdb,trpdb))
	if not os.path.exists(inpdb + '_a.out'):
                os.system(GEAR + '/enva_rec.v1.1 -a ' + trpdb + '.pdb > ' + trpdb + '_a.out')
        if not os.path.exists(inpdb + '_b.out'):
                os.system(GEAR + '/enva_rec.v1.1 -b ' + trpdb + '_het.pdb > ' + trpdb + '_b.out')
        if not os.path.exists(inpdb + '_m.out'):
                os.system(GEAR + '/enva_rec.v1.1 -m ' + trpdb + '.pdb ' + trpdb.split('.')[0] + 'B > ' + trpdb + '_m.out')
        if not os.path.exists('*.env'):
                os.system(GEAR + '/enva_rec.v1.1 -e ' + trpdb + '.pdb B')
#def tmp_ee(pdb):
#	trpdb = '_'.join(pdb.split('_')[:2])
#	if not os.path.exists('*.env'):
#		os.system(GEAR + '/enva_rec.v1.1 -e ' + trpdb + '_new.pdb B')

def red_rev(pdb):
	ppdb = pdb.split('.')#'.'.join(pdb.split('.')[])#'_'.join(pdb.split('_')[:2])
	pppdb = '.'.join([ppdb[0],ppdb[1],ppdb[3]])
	os.system('reduce -Trim ' + pdb + ' > ' + pppdb + '_red.pdb')
	os.system('csplit -f \'%s_red.pdb_\' %s_red.pdb \"/TER/\"'%(pppdb,pppdb))
	os.system('python /awork06-1/YKLee/py_script/pep_atom_revi.py ' + PDBLIB + '/' + inpdb + '_native.pdb ' + pppdb + '_red.pdb')
	os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py %s_rev.pdb A'%pppdb)
	os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py %s_rev.pdb B'%pppdb)
	os.system('cat %s_rev_A.pdb %s_rev_B.pdb > %s_rev_AB.pdb'%(pppdb,pppdb,pppdb))

def txt_writing(x):
        ff_list = []
        hh_list = []
        tt_acc = []
        tt_phi = []
        tt_psi = []
        tt_ac_dic = {}
        mat_lst = []
        lig_acc = []
        lig_num = []
        ave_acc = []
        tot_lig_acc = 0
        ave_lig_acc = 0
        tot_lig_num = 0
        ratio = 0
        zdx = 0
        pdbl = []
        tttt = []
        newl = []
	if not os.path.exists(x + '_rmsd.txt'):
		if x == 'traj_1':
			for i in sorted(glob.glob('*_AB_tr_B.pdb')):
                		os.system('/awork06-1/YKLee/c_script/rmsd_total ../../' + i.split('.')[0] + '_pep.pdb ' + i + ' bb >> ' + x + '_rmsd.txt') #PDBLIB + '/' + i.split('.')[0] + '_native.pdb ' + i + ' bb >> ' + x + '_rmsd.txt')
		elif x == 'v1':
			for i in sorted(glob.glob('*AB_tr_B.pdb')):
				os.system('/awork06-1/YKLee/c_script/rmsd_total ../../' + i.split('.')[0] + '_ref.pdb_rev_B.pdb ' + i + ' bb >> ' + x + '_rmsd.txt') 
        else:
		pass
	#rdf = pd.read_csv(x + '_rmsd.txt',sep='\s+',header=None)
        for f,g in zip(sorted(glob.glob('*a.out')),sorted(glob.glob('*m.out'))):
                sidx= []
                pdb = '_'.join([f.split('.')[0],str(x),f.split('.')[2].split('_')[0]])#'_'.join(f.split('.')[0].split('_')[:6])
                pdbl.append(pdb)
                adf = pd.read_csv(f,sep='\s+',skiprows=[0,0],header=None)
                mdf = pd.read_csv(g,sep='\s+',header=None)
                mdf_filter = mdf[mdf[17] == 1].iloc[:,9:16]
                for i in mdf_filter[9]:
                        tot_lig_acc = tot_lig_acc + float(i)
                        tot_lig_num += 1
                        if tot_lig_num > 0:
                                ave_lig_acc = tot_lig_acc/float(tot_lig_num)
                        else:
                                ave_lig_acc = 0
                ave_acc.append(ave_lig_acc)

                new_acc = {}
                num = 0
                for i in rfeats: #mdf_filter.T:
                        new_acc['AA_' + cc_dic[str(i)]] =sum([adf[adf[1] == int(i.split('_')[0])][adf[2] == i.split('_')[2]][adf[3] == i.split('_')[1]][11].values,1])
                for i in mdf_filter.T:
                        tidx = str(mdf_filter.T[i][13]) + '_' + str(mdf_filter.T[i][14]) + '_' + str(mdf_filter.T[i][15])
                        if tidx in rfeats:
                                tttt.append(tidx)
                                zdx = rfeats.index(tidx) +1
                                sidx.append(zdx)
                                num += 1
                newl.append(pd.DataFrame(new_acc))
                if not os.path.exists(pdb.split('.')[0] + '.ser'):
                        with open(pdb.split('.')[0] + '.ser','w') as W:
                                for sid in list(set(sidx)):
                                        W.write(str(sid)+'\n')
                os.system(GEAR +'/iskew ' + pdb.split('.')[0] + '.ser >> total_' + x + '_sk.txt')
                ratio = float(num)/float(len(rfeats))
                mat_lst.append(ratio)
        new_ac = pd.concat(newl)
	#new_ac.columns = rfeats1
	new_ac['total_rec_acc'] = new_ac.sum(axis=1)
        new_ac['ave_lig_acc'] = ave_acc
        new_ac['%Match'] = mat_lst
        #new_ac['rmsd'] = rdf[2]
        new_ac['PDB'] = pdbl
        nnew = new_ac.set_index('PDB').reset_index()

        nnew.to_csv(x + '_nac.txt',sep='\t',index=False)
        nn_list = []
        for f in sorted(glob.glob('*b.out')):
                with open(f,'r') as F:
                    n = 0
                    for line in F.readlines():
                        tt = line[7:13].strip() + '_' + line[17:21].strip() + '_' + line[13:17].strip()
                        if tt in rfeats:
                            n += 1
                    nn_list.append(n)
                with open(f,'r') as F:
                    hh_list.append(len(F.readlines()))
                ff_list.append(f.split('.')[0]+ '_' + str(x) + '_' + f.split('.')[2].split('_')[0])
        with open(x + '_hh.txt','w') as W:
                W.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
                for i,j,k in zip(ff_list,hh_list,nn_list):
                        W.write(str(i) + '\t' + str(j) + '\t' + str(k) + '\n')
	acc_columns = ['P%d'%i for i in range(1,seqlen+1)]
        phi_columns = ['PHI%d'%i for i in range(1,seqlen+1)]
        psi_columns = ['PSI%d'%i for i in range(1,seqlen+1)]
        pdbl = []
        for f in sorted(glob.glob('*.env')):
                pdb = '_'.join([sam.split('_')[0],str(x),f.split('.')[1].split('_')[0]])
                pdbl.append(pdb)
                with open(f,'r') as F:
                        acc = []
                        psi = []
                        phi = []
                        for line in F.readlines():
                                if line.find('chain') and line.startswith('ATOM'):
                                        envs = line.split()[9:]
                                        acc.append(envs[2])
                                        phi.append(envs[4])
                                        psi.append(envs[5])
                        tt_acc.append(pd.DataFrame(acc, index= acc_columns).T)
                        tt_phi.append(pd.DataFrame(phi, index = phi_columns).T)
                        tt_psi.append(pd.DataFrame(psi, index = psi_columns).T)
        pdbdf = pd.DataFrame(pdbl,columns = ['PDB'])
        acc_df = pd.concat(tt_acc,ignore_index=True).reindex()
        phi_df = pd.concat(tt_phi,ignore_index=True).reindex()
        psi_df = pd.concat(tt_psi,ignore_index=True).reindex()
        act_df = pd.concat([pdbdf,acc_df,phi_df,psi_df],axis=1)
        act_df.to_csv(x + '_ac_ct.txt',sep='\t',index=False)

def act(x,aa):
	#reduce and revise
	pdbs = sorted(glob.glob('*.pdb.*'))
	pool2 = multiprocessing.Pool(6)
	pool2.map(red_rev,pdbs)
	pool2.close()
	pool2.join()
	os.system('rm -rf *_01 *_00')
	os.system('rm -rf *.fasta')
	
	#sheba
	revpdbs = sorted(glob.glob('*rev_AB.pdb'))
	pool1 = multiprocessing.Pool(6)
	pool1.map(partial(sheba,b=aa),revpdbs)
	pool1.close()
	pool1.join()
	os.system('rm -rf *.trf')
	'''
	#reduce and revise
	trpdbs = sorted(glob.glob('*tr.pdb'))
	pool2 = multiprocessing.Pool(6)
	pool2.map(red_rev,trpdbs)
	pool2.close()
	pool2.join()
	os.system('rm -rf *tr.pdb')
	'''
	#enva
	trpdbs = sorted(glob.glob('*tr.pdb'))
	if x == 'traj_1' :
		pool3 = multiprocessing.Pool(6)
		pool3.map(enva_working,trpdbs)
		pool3.close()
		pool3.join()
	elif x == 'v1':
		pool3 = multiprocessing.Pool(6)
		pool3.map(original_enva,trpdbs)
		pool3.close()
		pool3.join()
	os.system('find ./ -name \"*_02\" -exec rm -rf {} \;')
	os.system('find ./ -name \"*_01\" -exec rm -rf {} \;')
	os.system('find ./ -name \"*_00\" -exec rm -rf {} \;')
	os.system('rm -rf *[0-9].pdb')
	'''
	#tmpdb = sorted(glob.glob('*new.pdb'))
	#tmpol = multiprocessing.Pool(5)
	#tmpol.map(tmp_ee,tmpdb)
	#tmpol.close()
	#tmpol.join()
	'''
	
	#produce txt file
	txt_writing(x)  #'traj_1')
	#os.system('find ./ -name \"*.out\" -exec rm {} \;')
	os.system('rm -rf *.ser')
	
	os.system('rm -rf *.fasta')
	#os.system('rm -rf *A.pdb')
	#os.system('rm -rf *B.pdb')
	#os.system('rm -rf *red.pdb')
	#os.system('rm -rf *rev.pdb')
	#os.system('rm -rf *.env')
	if x == 'traj_1' :
		os.system('ssh user1@10.1.5.5 \"mkdir -p /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf\"')
		os.system('tar cf - *het.pdb |ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf/\"')
		os.system('tar cf - *new.pdb |ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf/\"')
		os.system('tar cf - *a.out |ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf/\"')
		os.system('tar cf - *b.out |ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf/\"')
		os.system('tar cf - *m.out |ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf/\"')
		os.system('tar cf - *.env |ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/sim_conf/\"')

	else:
		os.system('find ./ -name \'*.pdb\' -exec rm {} \;')
	os.system('find ./ -name \'*.out\' -exec rm {} \;')
	os.system('find ./ -name \'*.pdb\' -exec rm {} \;')
	os.system('cp *.txt ../' + sam +'_energy_matrix/')
	os.chdir('../' + sam+'_energy_matrix/')
	
	df_rmsd = pd.read_csv(x + '_rmsd.txt',sep='\s+',header=None)
	df_rmsd.columns = ['pdb','reference','rmsd']
	df_hh = pd.read_csv(x + '_hh.txt',sep='\t')
	df_nac = pd.read_csv(x + '_nac.txt',sep='\t')
	df_sk = pd.read_csv('total_' + x +'_sk.txt',sep='\t',header=None).fillna('1.000000')
	df_sk.columns=['PDB','skewness','Class','Decision']
	total_df = [df_hh,df_nac,df_sk]
	df_final1 = reduce(lambda left,right: pd.merge(left,right,on=['PDB'],how='outer'),total_df)
	df_final1['rmsd'] = df_rmsd['rmsd']
	df_final2 = df_final1.set_index('rmsd').reset_index().set_index('PDB').reset_index()
	df_final2.to_csv(sam + '-' + x + '_total.txt',sep='\t',index=False,na_rep='-')
	os.chdir('../')
sam = sys.argv[1]
GEAR = '/awork06-1/neoscan_gear'
PDBLIB = '/awork06-1/YKLee/pdpdb/Neoscan_V2/Native/' + sam.split('_')[2] + '_native'
inpdb = sam.split('_')[0]
HLAclass = sam.split('_')[2]
seq = sam.split('_')[1]
seqlen = len(seq)
os.chdir(sam)
print(HLAclass)
print(seq)
#generate new reference PDB file in v1
os.system('csplit -f \'' + inpdb + '_ref.pdb_\' ' + inpdb + '_ref.pdb \'/TER/\'')
os.system('python /awork06-1/YKLee/py_script/pep_atom_revi.py ' + PDBLIB + '/' + inpdb + '_native.pdb ' + inpdb + '_ref.pdb')
os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py ' + inpdb + '_ref.pdb_rev.pdb B')
os.system('python /awork08/93_hong/rosetta_src_2019.14.60699_bundle/tools/protein_tools/scripts/clean_pdb.py ' + inpdb + '_ref.pdb_rev.pdb A')
os.system('cat ' + inpdb + '_ref.pdb_rev_A.pdb ' + inpdb + '_ref.pdb_rev_B.pdb > ' + inpdb + '_ref_new.pdb')
os.chdir('traj_1/')
try:
        if not os.path.exists(sam + '_energy_matrix'):
                os.mkdir(sam + '_energy_matrix')
except OSError:
        pass
rfeats = []
feats = pd.read_csv(PDBLIB + '/' + sam.split('_')[0] + '.out',sep='\s+',header=None)
feats_filter = feats[feats[17] == 1].iloc[:,13:16]
for i in feats_filter.values.tolist():
        rfeats.append(str(i[0]) + '_' + '_'.join(i[1:]))
rfeats1 = []
feats1_filter = feats[feats[17] == 1].iloc[:,12:16]
for i in feats1_filter.values.tolist():
	rfeats1.append(str(i[0]) + '_' + '_'.join(i[2:]))
cc_dic = {}
for i,j in zip(rfeats,rfeats1):
	cc_dic[i] = j

#activation part
os.chdir('pdb_from_prod/')
os.system('cp ../../feat.list .')
act('traj_1',PDBLIB + '/' + inpdb + '_native.pdb')
print(HLAclass)
print(seq)
os.system('tar cf - ' + sam + '_energy_matrix | ssh 10.1.5.5 \"tar xpfm - -C /awork05/d01/Neo_out/v2/' + str(seqlen) + 'mer/' + HLAclass + '/' + inpdb + '_final_res/\"')
#original calculation parts
#os.chdir('pdb_from_prod/')
#act('v1','/home/user1/' + sam + '/' + inpdb + '_ref_new.pdb')
