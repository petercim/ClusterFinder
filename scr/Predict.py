###################################################################################################################################
#                                                                                                                                 #
#                                                       ClusterFinder Library                                                     #
#                                                         Peter Cimermancic                                                       #
#                                                               2010                                                              #
#                                                                                                                                 #
###################################################################################################################################



import sys
from math import log
from random import random
import numpy
from HMMs import *
import pickle, copy
from operator import itemgetter


def do_HMM(observations,A,B,index_dict,start_p,graph=False,report=False,outs='',method='HMM',path='/pico1/home/peterc/Documents/Fischbach_Lab/',name='bla'):

	'''
	This is an interface between HMM module and users data. It takes in Python array of Pfam IDs ordered by their genomic location, and returns an array
	of gene cluster probabilities. Optionally, it can also write a text report and plot a lousy probabilites versus pfam position graph. 

	Inputs:
		Mandatory:
		- observations: Python array of Pfam IDs (example: ['PF00001','PF00550','PF00501',...])
		- A: Numpy array of transition probabilities (read HMM module documentation for more info)
		- B: Numpy array of emission probabilites (read HMM module documentation for more info)
		- index_dictio: dictionary of Pfam indeces for positions in A
		- start_p: Numpy vector of start probabilities
		Optional:
		- graph(0/1): plot graph (0 by default)
		- report(0/1): write report (0 by default)
		- outs: Pyton array of arrays with information about each Pfam you would like to write in output,
			    a probability will be append to each of the sub-arrays, tab-delimited, and written in the output file
		- method: HMM for Forward-Backward algorithm, viterbi for Viterbi alogirthm (HMM by default)
		- path: path to input and output files
		- name of the output files (graph and/or report)
	'''

	index_dictio = copy.deepcopy(index_dict)
	    
	reals = None
	if method=='viterbi':
		print '\tPredicting Clusters...'
		PREDICTED = list((numpy.array(viterbi(A,B,start_p,observations,index_dictio))-1)*(-1))
		print '\tCluster predicted!\n'

	elif method == 'HMM':   
		print '\tPredicting Clusters...'
		PREDICTED = list(forward_backward(A,B,start_p,observations,index_dict,method='BF'))
		print '\tCluster predicted!\n'
	
	else:
		print "\nMethod: '%s' unknown. Please choose HMM or viterbi\n" % method
		sys.exit(1) 


	if graph == True and method=='HMM':
		ploter(PREDICTED,750,name,path)  # output name as an argument
	if graph == True and method == 'viterbi':
		ploter(PREDICTED,750,name,path,real_ranges=reals)  # output name as an argument
		# in cases, they won't be real any more but rather predicition of viterbi alg.

	if report == True:
		output = open(path+name+'.out','w')
		if len(outs) == len(PREDICTED):
			for n,i in enumerate(PREDICTED):
				outs[n].append(PREDICTED[n])
				lines = '\t'.join([str(i) for i in outs[n]]) + '\n'
				output.write(lines)
		else:
			print 'OutputWriteError: PREDICTED array length does not match wigh Genome class!'
			sys.exit(1)
		output.close()

	print '%s Done!' % name
	return PREDICTED

def ploter(X,n,name,path,real_ranges = None):
	'''
	Graph plotting module. Very basic, so feel free to modify it.
	'''
	try: import matplotlib.pylab as plt
	except ImportError: 
		print '\nMatplotlib is not installed! Either install it, or deselect graph option.\n'	
		return 1

	lll = 1+len(X)/n
	fig = plt.figure(figsize=(16,14),dpi=100)
	for x in xrange(0,len(X),n):
		iii = 1+x/n
		ax = fig.add_subplot(lll,1,iii)
       
		if real_ranges != None:
			for p in real_ranges:
				if p[0] in range(x,x+n) and p[1] in range(x,x+n): ax.axvspan(p[0]%(n), p[1]%(n), facecolor='g', alpha=0.5)
				elif p[0] in range(x,x+n) and p[1] > x+n: ax.axvspan(p[0]%(n), n, facecolor='g', alpha=0.5)
				elif p[0] < x and p[1] in range(x,x+n): ax.axvspan(0, p[1]%(n), facecolor='g', alpha=0.5)
				elif p[0] < x and p[1] > x+n: ax.axvspan(0, n, facecolor='g', alpha=0.5)

		ax.plot(X[x:x+n],'r-')

		ax.set_xlim(0.,n)
		ax.set_ylim(0.,1.)
		ax.set_xticklabels(range(x,x+750,100)) #(x,x+n/5,x+2*n/5,x+3*n/5,x+4*n/5,x+5*n/5) )
                
	plt.savefig(path+'HMM_'+name+'.png')
	plt.close()


def repeat_filter(observations):
	'''
	Takes Python array of Pfam IDs in ant filters it for repeats.
	'''
	repeats = ['PF07721', 'PF05593', 'PF07719', 'PF00515', 'PF00132', 'PF03130', 'PF01839', 'PF01816', 'PF07720', 'PF00400', 'PF05594', 'PF07661', 'PF02985', 'PF06049', 'PF08238', 'PF06696', 'PF00353', 'PF02412', 'PF00023', 'PF02071', 'PF03991',   	'PF01469', 'PF07676', 'PF00514', 'PF00904', 'PF07634', 'PF02370', 'PF03335', 'PF01851', 'PF04728', 'PF06715', 'PF03373', 'PF04680', 'PF00805', 'PF04508', 'PF07918', 'PF01535', 'PF01011', 'PF05017', 'PF06671', 'PF00818', 'PF03406', 'PF00399',   	'PF09373', 'PF01744', 'PF01436', 'PF01239', 'PF05906', 'PF03729', 'PF00404', 'PF04022', 'PF02363', 'PF02524', 'PF07981', 'PF02095', 'PF00414', 'PF00560', 'PF05001', 'PF02162', 'PF01473', 'PF05465', 'PF02493', 'PF03578', 'PF08043', 'PF06392',   	'PF07142', 'PF08309', 'PF02184']

	return [o for o in observations if o not in repeats]

# TEST
#print do_HMM(repeat_filter(['PF07719','PF00001','PF00550','PF00501']),A,B,index_dict,start_p,method='HMM',graph=0,name='testtest')

