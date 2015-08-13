# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 17:58:07 2015

@author: andrea.insabato@upf.edu

Decision-making attractor neural network with integrate & fire neurons with AMPA, GABA and NMDA synapses
from Wang (2002).

TODO: check consistency with Wang data
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

unitsCheck = False  
if not unitsCheck: import brian_no_units
from brian import *

clear()
#from brian.globalprefs import *
#set_global_preferences(useweave=True)

import numpy
import time 
import gzip
import cPickle as pickle
import random as rnd
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import argparse

parser = argparse.ArgumentParser(description="run simulations of a decision-making network composed of rate models")
parser.add_argument("-W", "--adjacency_matrix_file", help="the file where the connectivity matrix is stored (the file must be in .pickle.gz format)")
args = parser.parse_args()

# random seed
#seed = int(time.time())
# fixed seed: 1423208210
seed = 1423208210
numpy.random.seed(seed)

exec(open('defs.py', 'rb').read()) # Constants, etc.

#---- function to create connectivity matrix ------
def create_W(pre_nrns, post_nrns, sparse_fact):
    # hubs_list is a list with the indices of units that will be hubs
    # hubness is a list indicating the proportion of links for each hub (do not confuse it with the hubness measure in GAlib!)
    W = numpy.ones([pre_nrns, post_nrns])
    # set some connections to 0
    #  eliminate 1-sparse_fact links in the whole network (resulting in gaussian degree distr)
    tmp = numpy.reshape(W, -1)
    tmp[rnd.sample(range(size(W)), int(size(W)*(1-sparse_fact)))] = 0
    W = numpy.reshape(tmp, shape(W))
    return W
#--------------------------------------------------


def run_single_trial(trial_num):
    
    #---- net op to calculate the NMDA input ----------
    @network_operation(when='start') 
    def sNMDA_contribution():
        s_NMDA_e = numpy.dot(Pe.s_nmda_out, W_ee)
        s_NMDA_i = numpy.dot(Pe.s_nmda_out, W_ei)    
        # NMDA between exc groups
        Pe.s_nmda_in = s_NMDA_e
        # NMDA from exc to inh group
        Pi.s_nmda_in = s_NMDA_i
    #--------------------------------------------------

    #---- net op to stop simulation when a decision is met ----------
    @network_operation() 
    def stops_ifDecided():
        if defaultclock.t > 400.*ms:
            rateL = rate_monitors[1].smooth_rate(width=50*ms,filter='flat')[-20]
            rateR = rate_monitors[2].smooth_rate(width=50*ms,filter='flat')[-20]
            X = abs( (rateL - rateR) / (rateL + rateR) )
            if X > 0.75:
                stop()
    #----------------------------------------------------------------

    defaultclock.t = 0.0*ms  # start the time at 0 in order to compare spike times
    defaultclock.dt = 0.02*ms
    duration = 1.*second
    
    dt = defaultclock.dt
    
    #----------- neuron's equations -------------------
    eqs = ''' 
    dv/dt = (-gm*(v-VL) - \
             (g_ampa_ext*(v-VE)*s_ampa_ext + \
              g_ampa_rec*(v-VE)*s_ampa_rec + \
              g_nmda*(v-VE)*s_nmda_in/(1+exp(-0.062*v/(1*mV))/3.57) + \
              g_gaba*(v-VI)*s_gaba) \
            )/Cm: mvolt  
    ds_ampa_ext/dt = -s_ampa_ext/(t_ampa) : 1
    ds_ampa_rec/dt = -s_ampa_rec/(t_ampa) : 1
    s_nmda_in : 1
    ds_nmda_out/dt = -s_nmda_out/t_nmda_decay + alpha_nmda*x*(1-s_nmda_out) : 1
    ds_gaba/dt = -s_gaba/(t_gaba) : 1
    dx/dt = -x/(t_nmda_rise) : 1
    gm : siemens
    g_ampa_ext : siemens
    g_ampa_rec : siemens
    g_nmda : siemens
    g_gaba : siemens
    Cm: farad
    '''
    #--------------------------------------------------
    
    print 'Populations...' 
    # All neurons
    P = NeuronGroup(1000, eqs, threshold=Vthr, reset=Vreset, 
                    refractory=r_[numpy.array([tref_e]).repeat(800),
                                  numpy.array([tref_i]).repeat(200)],
                    unit_checking=unitsCheck,
                    order=2) # RK2 is closer to Heun
                    
    P.v = Vreset+numpy.random.random(1000)*(Vthr-Vreset);
    P.s_ampa_ext, P.s_ampa_rec, P.s_gaba, P.s_nmda_in, P.s_nmda_out, P.x = \
            0, 0, 0, 0, 0, 0
    
    # Excitatory and inhibitory groups
    Pe = P.subgroup(800)
    Pi = P.subgroup(200)
    
    # Constant state variables
    Pe.gm, Pe.g_ampa_ext, Pe.g_ampa_rec, Pe.g_nmda, Pe.g_gaba, Pe.Cm = \
            gm_e, g_ampa_ext_e, g_ampa_rec_e, g_nmda_e, g_gaba_e, Cm_e
    Pi.gm, Pi.g_ampa_ext, Pi.g_ampa_rec, Pi.g_nmda, Pi.g_gaba, Pi.Cm = \
            gm_i, g_ampa_ext_i, g_ampa_rec_i, g_nmda_i, g_gaba_i, Cm_i
    
    Pl = Pe.subgroup(120)
    Pr = Pe.subgroup(120)
    Pns = Pe.subgroup(560)
        
    Pp = PoissonGroup(1000, rates = 2.4*kHz)
    Ppl = PoissonGroup(120, rates = lambda t: (t<.5*second and 0) or 50)
    Ppr = PoissonGroup(120, rates = lambda t: (t<.5*second and 0) or 50)
    
    print 'Connections...'
    struct_type = 'sparse'
    # external synapses
    Cp = IdentityConnection(Pp, P,'s_ampa_ext', weight=1.0)
    Cpl = IdentityConnection(Ppl, Pl,'s_ampa_ext', weight=1.0)
    Cpr = IdentityConnection(Ppl, Pr,'s_ampa_ext', weight=1.0)
    # take connectivity matrix from file if specified or create one
    if args.adjacency_matrix_file:
    	print "using connectivity matrix in {}".format(args.adjacency_matrix_file)
    	f = gzip.open(args.adjacency_matrix_file, 'r')
    	W = pickle.load(f) 
    	f.close()
    else:
	print "WARNING: connectivity file not specified! Creating a random graph with density {}".format(sparse_fact)
        # E-->E connection matrix
        W_ee = create_W(800, 800, 1.)
        W_ss = create_W(240, 240, sparse_fact)
        W_ss[0:120,0:120] *= wp_sparse
        W_ss[120:,120:] *= wp_sparse
        W_ss[0:120,120:] *= wm_sparse
        W_ss[120:,0:120] *= wm_sparse
        W_ee[0:240,0:240] = W_ss
        W_ee[240:,0:240] *= wm
        # E-->I connection matrix
        W_ei = create_W(800, 200, 1.)
        # I-->I connection matrix
        W_ii = create_W(200, 200, 1.)
        # I-->E connection matrix
        W_ie = create_W(200, 800, 1.)
    
    # E-->E
    Cee = Connection(Pe, Pe, 's_ampa_rec', weight = W_ee,
                     structure=struct_type, sparseness = 1)
    # E-->I
    Cei = Connection(Pe, Pi, 's_ampa_rec', weight = W_ei,
                     structure=struct_type, sparseness = 1)
    # I-->I
    Ciis = Connection(Pi, Pi, 's_gaba', weight=W_ii, 
                       structure=struct_type, sparseness= 1)
    # I-->E
    Cies = Connection(Pi, Pe, 's_gaba', weight=W_ie, 
                       structure=struct_type, sparseness= 1)
    
    # this fictitious connection increases presynaptic nmda x variable by 1.0
    # when a neuron emits a spike. The alternative is to modify the reset procedure
    # but this, in turn, doesn't allow to set different refractory periods for
    # each neuron, so NeuronGroup P need to be separated in exc and inh. Execution
    # time for the two alternatives is practically the same.
    C_nmda = IdentityConnection(Pe, Pe, 'x', weight = 1.0)
    
    print 'Monitors...'
    
    R = SpikeMonitor(P)
    rate_monitors = [PopulationRateMonitor(pool, 1*ms) for pool in [Pi, Pl, Pr, Pns]]
    
    print 'Running... trial ' + str(trial_num)
    
    run(duration, report='text', report_period = 30*second)
    
    figure()
    for i in range(len(rate_monitors)):
        plot(rate_monitors[i].times, rate_monitors[i].smooth_rate(width=50*ms,filter='flat'))
    legend(('Inh','L','R','Ns'))
    show()
    return R

dataman = DataManager('prova_ANN')
run_tasks(dataman, run_single_trial, range(1))
S = dataman.values()
fname = 'prova.pickle.gz'
with gzip.open(fname,'wb') as f:
    pickle.dump(S, f)

