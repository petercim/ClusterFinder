###################################################################################################################################
#                                                                                                                                 #
#                                                   HIDDEN MARKOV MODEL ALGORITHMS                                                #
#                                                         Peter Cimermancic                                                       #
#                                                  Reference: Numerical Recepies in C                                             #
#                                                               2010                                                              #
#                                                                                                                                 #
###################################################################################################################################


import numpy
import copy

BIG = 1.e20
BIGI = 1./BIG
STEP = 0.

def forward_backward(A,B,start_p,observations,DICT,method='BF'):
    '''
    mstat, nobs, ksym ... numbers of states obervations, and symbols (different domains)
    alpha, beta, pstate ... matrices:
        The observable symbol probability
        distributions are represented by an N X M matrix where M is the number of observation
        symbols.  

                  |a_11 a_12 ... a_1N|             |b_11 b_12 ... b_1M|            |s1|
                  |a_21 a_22 ... a_2N|             |b_21 b_22 ... b_2M|            |s2| 
              A = | .    .        .  |         B = | .    .        .  |   start_p =| .| 
                  | .         .   .  |             | .         .   .  |            |  |
                  |a_N1 a_N2 ... a_NN|             |b_N1 b_N2 ... b_NM|            |sN|
          
               a_ij = P(q_t = S_j|q_t-1 = S_i)       b_ik = P(v_k at t|q_t = S_i)
           where q_t is state at time t and v_k is k_th symbol of observation sequence.

    Dishonest casino example:

           A = numpy.array([[F->F, F->L],     B = numpy.array([[F1, F2, F3, F4, F5, F6],    start_p =numpy.array([[0.6, 0.4]])
                            [L->F, L->L]])                     [L1, L2, L3, L4, L5, L6]])             
    DICT is a dictionary of alphabates and B array indexes.

    '''
    DICT_copy = copy.deepcopy(DICT)
    global BIG,BIGI
    mstat, nobs, ksym = len(A), len(observations),len(B[0,:]) 
    #obervations = numpy.array([1,2,3,2,4,6,2,])
    alpha = numpy.zeros((nobs+1,mstat),dtype=numpy.float64())

               ### FORWARD ALGORITHM ###
    for i in xrange(mstat):
        try:
            alpha[0][i] = B[i][DICT[observations[0]]]
        except KeyError:
            alpha[0][i] = 1.
            DICT_copy[observations[0]] = ksym
            B = numpy.hstack((B,[[1./ksym],[1./ksym]]))
            ksym += 1
    arnrm = [0]  
    for t in xrange(1,nobs):
        asum = 0.
        for j in xrange(mstat):
            sum = 0
            for i in xrange(mstat):
                try:
                    sum += alpha[t-1][i]*A[i][j]*B[j][DICT[observations[t]]]  
                except KeyError: 
                    '''
                    Since PFAM domain doesn't exist in neither DICT or B, we will create room for it.
                    Since we don't have good aproximation of what its frequency should be, we will
                    ignore it at calculation of <sum>, but the frequency for better guess at another cycle
                    of Baum-Welch is set to frequency of a single domain set of all of the domains.
                    '''
                    sum += alpha[t-1][i]*A[i][j]*1.
                    DICT_copy[observations[t]] = ksym
                    B = numpy.hstack((B,[[1./ksym],[1./ksym]]))
                    ksym += 1
            alpha[t][j] = sum
            asum += sum
        
        arnrm.append(arnrm[-1])
        if asum < BIGI:                # normalization to avoid float underflow
            arnrm[-1] += 1
            for j in xrange(mstat):
                alpha[t][j] *= BIG

    arnrm = numpy.array(arnrm)
    #alpha = alpha[1:,:] # cut off start_p vector

   
                   ### BACKWARD ALGORITHM ###
    beta = numpy.zeros((nobs,mstat),dtype=numpy.float64())   # initialization
    beta[-1] = numpy.ones((1,2),dtype=numpy.float64())
    brnrm = numpy.zeros(nobs,dtype=numpy.int32())        # begin passages
    for t in xrange(nobs-2,-1,-1):
        bsum = 0.
        for i in xrange(mstat):
            sum = 0.
            for j in xrange(mstat):
                try:
                    sum += A[i][j]*B[j][DICT[observations[t+1]]]*beta[t+1][j]
                except KeyError: sum += A[i][j]*1.*beta[t+1][j]
            beta[t][i] = sum
            bsum += sum
        
        brnrm[t] = brnrm[t+1]
        if bsum < BIGI:
            brnrm[t] += 1
            for j in xrange(mstat):
                beta[t][j] *= BIG
                
                ### CALCULATE LIKELIHOOD ###
    lhood = 0.                            # overall likelihood
    for i in xrange(mstat):
        lhood += alpha[0][i]*beta[0][i]
    lrnrm = arnrm[0] + brnrm[0]
    while lhood < BIGI:                  # again, to avoid overflow problem
        lhood *= BIG
        lrnrm += 1
    
    pstate = numpy.zeros((nobs,mstat),dtype=numpy.float64())
    for t in xrange(nobs):
        sum = 0.
        for i in xrange(mstat):
            pstate[t][i] = alpha[t][i]*beta[t][i]
            sum += pstate[t][i]
        for i in xrange(mstat):
            pstate[t][i] /= sum
    
    if method == 'BF':        
        return pstate[:,0]
    elif method == 'BW':
        return baum_welch(A,B,start_p,alpha,beta,arnrm,brnrm,lrnrm,lhood,observations,DICT_copy)


def baum_welch(A,B,start_p,alpha,beta,arnrm,brnrm,lrnrm,lhood,observations,DICT):
    '''
    Baum-Welch procedure
    '''
    global STEP,BIGI
    
    olda = copy.deepcopy(A)
    oldb = copy.deepcopy(B)
    oldstart_p = copy.deepcopy(start_p)
    mstat, nobs, ksym = len(A), len(observations),len(B[0,:]) 
    bnew = numpy.zeros((mstat,ksym),dtype=numpy.float64())
    
    powtab = numpy.ones((10),dtype=numpy.float64())
    for i in xrange(len(powtab)):
        powtab[i] = numpy.power(BIGI,i-6)
    for i in xrange(mstat):
        denom = 0.
        for t in xrange(nobs):
            term = (alpha[t][i]*beta[t][i]/lhood)*powtab[arnrm[t] + brnrm[t] - lrnrm + 6]
            denom += term
            bnew[i][DICT[observations[t]]] += term
        for j in xrange(mstat):
            num = 0.
            for t in xrange(nobs-1):
                num += alpha[t][i]*B[j][DICT[observations[t+1]]]*beta[t+1][j]*powtab[arnrm[t] + brnrm[t+1] - lrnrm + 6]/lhood
            A[i][j] *= (num/denom)
        for k in xrange(ksym):
            bnew[i][k] /= denom
        start_p[i] = denom/nobs   # new way to calculate prior probabilities
    B = bnew
    
    # calculate mean of changes in A,B, and start_p together
    meana = numpy.sum(numpy.abs(A-olda)) / numpy.size(A)
    meanb = numpy.sum(numpy.abs(B-oldb)) / numpy.size(B) 
    meanstart_p = numpy.sum(numpy.abs(start_p-oldstart_p)) / numpy.size(start_p)
    mean = (meana + meanb + meanstart_p) / 3
    STEP += 1
    #print STEP,mean
    #print A
    #print B
    #print start_p
    #print
    #
    # termination test: if steps >= x or deltas <= y: return pstate
    print 'STEP: %i mean: %.9f' % (STEP,mean)
    if mean < 1e-6 or STEP == 100:
        return forward_backward(A,B,start_p,observations,DICT,method='BF')
    else:
        return forward_backward(A,B,start_p,observations,DICT,method='BW')
        
    

def viterbi(A,B,start_p,observations,DICT):
    '''
    Calculate the most probable path, and return it!
    WARNING: Not well tested!


    mstat, nobs, ksym ... numbers of states obervations, and symbols (different domains)
    alpha, beta, pstate ... matrices:
        The observable symbol probability
        distributions are represented by an N X M matrix where M is the number of observation
        symbols.  

                  |a_11 a_12 ... a_1N|             |b_11 b_12 ... b_1M|            |s1|
                  |a_21 a_22 ... a_2N|             |b_21 b_22 ... b_2M|            |s2| 
              A = | .    .        .  |         B = | .    .        .  |   start_p =| .| 
                  | .         .   .  |             | .         .   .  |            |  |
                  |a_N1 a_N2 ... a_NN|             |b_N1 b_N2 ... b_NM|            |sN|
          
               a_ij = P(q_t = S_j|q_t-1 = S_i)       b_ik = P(v_k at t|q_t = S_i)
           where q_t is state at time t and v_k is k_th symbol of observation sequence.

    Dishonest casino example:

           A = numpy.array([[F->F, F->L],     B = numpy.array([[F1, F2, F3, F4, F5, F6],    start_p =numpy.array([0.6, 0.4])
                            [L->F, L->L]])                     [L1, L2, L3, L4, L5, L6]])             
    DICT is a dictionary of alphabates and B array indexes.
    '''

    mstat, nobs, ksym = len(A), len(observations),len(B[0,:]) 
    T = {}
    
    for m in xrange(mstat):
        T[m] = (numpy.log(start_p[m]),[m],start_p[m])
    for t in xrange(0,nobs):
        U = {}
        for j in xrange(mstat):
            total = 0
            valmax = None
            argmax = None
            for i in xrange(mstat):
                (prob,v_path,v_prob) = T[i]
                try:
                    p = numpy.log(A[i][j]*B[j][DICT[observations[t]]])
                except KeyError:
                    p = numpy.log(A[i][j]*1.)
                prob += p
                v_prob += p
                total += prob            
                if v_prob > valmax:
                    argmax = v_path + [j]
                    valmax = v_prob
            U[j] = (total,argmax,valmax)
        T =  U
      
    ## apply sum/max to the final states:
    total = 0.
    argmax = None
    valmax = None
    for t in xrange(mstat):
        (prob, v_path, v_prob) = T[t]
        total += prob
        if v_prob > valmax:
            argmax = v_path
            valmax = v_prob
    return argmax[:-1]



### TEST ###
'''
A = numpy.array([[0.80,0.20],
                 [0.4,0.6]])
B = numpy.array([[1., 1./6, 1./6, 1./6, 1./6, 1./6],
                 [1./6, 1./6, 1./6, 1./6, 1./6, 1./6]])
start_p =numpy.array([0.6, 0.4])
observations = '4462651645252332266216166466424135546113323166315362545435552656455212616666556351512621535451641666654642222124216615544665646232441432445533132411161246335115346155254616552142154435451421366666246424411264236663656336463354634651166254254652564535625452635336445653153636331363462666652611656666646211266261646663662136614344366142666412244613165513161642464226636615155331521643415236616646364642355626114141536514423232354653356243614166634163154524566453463665614454322525143215313241113361464452366166651265434426666666661165354665351666166456341666114644345626546362266265365441624453424524315226456265334442436262212661345112221266564354666626311664166664546511166653261665126644566646364661621635355464612634165254444662265261666236646161232451234234256136511222666166423266245146445461636511551642254612156635166655433253145645463522332454266366662442521142244625341455513151365666431343666661664435165665666156666451333331566266546442311665655431134462514226255451252514352144432145622636461443134223244441145123625326313636135634666614211143165444613446566636146266644563365611335546665232621245666421332246156614561552255314566415356421541314663116416656316626566164666412241415636264236122134541356323463232131655624641554433366666161656611566512142123642666415533234345133613154513246143236552522145431541464645445461466266166662626542665154425524626652622255244124526425225154135441613431261236331562425646516645242313345452265546665123313665466446263355324666513665135353441234363644551653255464245515234633411311261416211326365644615645345526346614342645141411522333134141136161362166366324234646422416456555163323655556462115534465423354513462623355243342466611134245666162411413421222412224155651444226421663156511461666113521253113566663253344156566312663561666156665512211355241116661566662566666662126511465323116666112236563633435142612346251143421535156125642311132532461541413553644655462441611453222311514561426535256665446113114546121145112462334131235541132545436321442652'
REAL = 'LLLLLLLLLLLLLFLLFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLLLLFFFFFFFFFFFFFFFFFFFFFLLLLLL'
observation = numpy.array([i for i in observations])
DICT = {'1':0,'2':1,'3':2,'4':3,'5':4,'6':5}

reals = []
for i in xrange(len(REAL)-1):
    
    if (REAL[i] == 'F' and REAL[i-1] != 'F') or (REAL[i] == 'F' and i==0):
        sten = [0,0]
        sten[0] = i
    elif REAL[i] != 'F' and REAL[i-1] == 'F': 
        sten[1] = i
        reals.append(sten)
#PRED = forward_backward(A,B,start_p,observations,DICT,method='BW')
PRED = viterbi(A,B,start_p,observations,DICT)
print PRED
'''
