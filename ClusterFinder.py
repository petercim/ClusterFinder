##############################################################################
#                    ClusterFinder - main script                             #
#                       Peter Cimermancic                                    #
##############################################################################


PATH = './scr'
import sys
sys.path.append(PATH)
from Predict import *
from FindClusters import *
import re


organism_name = 'example_org'
input_file = 'example_input.txt'

# --- read in the HMM model data
PATH_TO_FREQUENCIES = './freqdata/' # link to where the frequencies are
PATH_TO_DATA = './data/'

out4 = open(PATH_TO_FREQUENCIES + 'TP_arr_A_latest.pkl','rb')
out5 = open(PATH_TO_FREQUENCIES + 'NewTS_all_B_reduced_6filter.pkl','rb')
out6 = open(PATH_TO_FREQUENCIES + 'SP_arr.pkl','rb')
out7 = open(PATH_TO_FREQUENCIES + 'NewTS_all_B_index.pkl','rb')

A=pickle.load(out4)
B=pickle.load(out5)
start_p=pickle.load(out6)
index_dict=pickle.load(out7)

out4.close()
out5.close()
out6.close()
out7.close()


# --- read in an input file
'''
Costumize according to your input.
'''

data = open(PATH_TO_DATA+input_file)
D = data.readlines()
data.close()

DATA = []
for d in D:
	d = d.strip().split('\t')
        
	gene,pfamID = d[5],d[13].replace('pfam','PF')
	genstart,genestop = int(d[6]),int(d[7])
	try: pfamstart,pfamstop = int(d[11]),int(d[12])
	except ValueError: pfamstart,pfamstop = genstart,genestop
	chain = d[3]
	
	# - if gene has a Pfam annotation
	if pfamID != 'n/a': DATA.append([ chain,gene,genstart,genestop,pfamstart,pfamstop,pfamID ])
	# - else
	else: DATA.append( [ chain,gene,genstart,genestop,pfamstart,pfamstop,'Q' ] )


# --- sort Pfams as they appear in the genome
PfamSorted = sorted(DATA,key = itemgetter(0,2,4))

# --- filter out the repeats (optional)
repeats = set(['PF07721', 'PF05593', 'PF07719', 'PF00515', 'PF00132', 'PF03130', 'PF01839', 'PF01816', 'PF07720', 'PF00400', 'PF05594', 'PF07661', 'PF02985', 'PF06049', 'PF08238', 'PF06696', 'PF00353', 'PF02412', 'PF00023', 'PF02071', 'PF03991', 'PF01469', 'PF07676', 'PF00514', 'PF00904', 'PF07634', 'PF02370', 'PF03335', 'PF01851', 'PF04728', 'PF06715', 'PF03373', 'PF04680', 'PF00805', 'PF04508', 'PF07918', 'PF01535', 'PF01011', 'PF05017', 'PF06671', 'PF00818', 'PF03406', 'PF00399', 'PF09373', 'PF01744', 'PF01436', 'PF01239', 'PF05906', 'PF03729', 'PF00404', 'PF04022', 'PF02363', 'PF02524', 'PF07981', 'PF02095', 'PF00414', 'PF00560', 'PF05001', 'PF02162', 'PF01473', 'PF05465', 'PF02493', 'PF03578', 'PF08043', 'PF06392', 'PF07142', 'PF08309', 'PF02184'])

PfamSorted = [i for i in PfamSorted if i[-1] not in repeats]


# - check that there are Pfam domain in the data
pfpatt = re.compile(r"""PF+\d+\d+\d+\d+\d""")
if len([re.match(pfpatt,p[-1]) for p in PfamSorted])==0:
	print 'There is no Pfam domain in you dataset!'
	exit()

# --- run the HMM
do_HMM([i[-1] for i in PfamSorted],A,B,index_dict,start_p,graph=False,report=True,outs=PfamSorted,\
	method='HMM',path=PATH_TO_DATA,name=organism_name)


# --- filter for biosynthetic domains
data = open(PATH_TO_FREQUENCIES+'biosynthetic_pfams.txt','r')
BioPfams = set([i.strip().split('\t')[0] for i in data.readlines()])
data.close()


# --- run clustering method
Clusters, OrgName = get_clusters(PATH_TO_DATA+organism_name+'.out',X=2,TRSH=0.2)
CLUSTERS = [] # store the new clusters
for clstr in Clusters:
	pfams =  set([i[-4] for i in clstr])
	if len(pfams - BioPfams)>=3: # - requires at least three biosynthetic domains
		CLUSTERS.append(clstr)


# --- write the data (in ClusterFinder format)
output = open(PATH_TO_DATA+organism_name+'.clusters.out','w')
for x,clstr in enumerate(CLUSTERS):
	for domain in clstr:
		output.write('\t'.join([str(j) for j in [organism_name+'_'+str(x)]+domain])+'\n')
output.close()



