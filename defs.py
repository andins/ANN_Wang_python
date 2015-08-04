# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 17:01:27 2015

@author: andrea insabato
"""

# === Parameters determining model dynamics =====

sparse_fact = .2

# pyramidal cells 
Cm_e        = 0.5*nF    # total capacitance
gm_e        = 25.0*nS   # total leak conductance 
tref_e      = 2.0*ms    # refractory time 
 
# interneuron cells 
Cm_i        = 0.2*nF    # total capacitance 
gm_i        = 20.0*nS   # total leak conductance 
tref_i      = 1.0*ms    # refractory time

# Common values
Vthr      = -50.0*mV  # threshold potential
Vreset    = -55.0*mV  # reset potential 
VL        = -70.0*mV  # leak reversal potential 

VE          = 0.0*mV    # synaptic reversal potential for excitation
VI          = -70.0*mV  # synaptic reversal potential for inhibition

# AMPA receptor (AMPAR) 
t_ampa      = 2.0*ms    # exponential decay time constant  
g_ampa_ext_e= 2.1*nS  # maximum conductance from external to pyramidal cells 
g_ampa_ext_i= 1.62*nS  # maximum conductance from external to interneurons
g_ampa_rec_e= 0.1*nS
g_ampa_rec_i= 0.08*nS

# NMDA receptor (NMDAR) 
t_nmda_rise = 2.0*ms    # controls the rise time of NMDAR channels 
t_nmda_decay= 100.0*ms  # decay time of NMDA currents 
g_nmda_e    = 0.33*nS
g_nmda_i    = 0.26*nS
alpha_nmda  = 0.5*kHz   # controls saturation properties of NMDAR channels
#beta_nmda = 0.062/mV    # hardcoded (Brunel&Wang's formulation)
#gamma_nmda= 0.28011     # hardcoded (this values are in the neuron's eqs)

# GABA receptor (GABAR) 
t_gaba      = 5.0*ms   # exponential decay time constant
g_gaba_e    = 2.6*nS
g_gaba_i    = 2.*nS
 
 
## variables
#
NI = 200
NE = 800
NL = 120
NR = 120
Nns = 560

wp = 1.7
wm = ( (1-.15*wp)/(1-.15) ) 
wp_sparse = wp / sparse_fact
wm_sparse = wm / sparse_fact
