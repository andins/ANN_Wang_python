""" TODO: check for the existence of clusters: how well connected is the network? """

import numpy as np
import random as rnd
from pylab import *
import gzip
import cPickle as pickle
#import gamodels
import galib
import argparse

parser = argparse.ArgumentParser(description="create a connectivity matrix for ANN_wang2002")
parser.add_argument("-o", "--out_file",
        help="the file name where the graph will be saved (extention .pickle.gz will be added)",
        default="W")
parser.add_argument("-n", "--num_neurons", help="number of neurons in the network",
        type=int, default=240)
parser.add_argument("-d", "--density", help="the density of the network",
        type=float, default=.2)
parser.add_argument("-H", "--hubs_list", help="a list of neurons that will have more links",
        type=int, nargs='+', default=[])
parser.add_argument("-D", "--hubs_density", help="a list of densities for each hub in the network",
        type=float, nargs='+', default=[])
args = parser.parse_args()

def create_W(num_nrns, density, hubs_list, hubs_density, seed='none'):
    # hubs_list is a list with the indices of units that will be hubs
    # hubs_density is a list indicating the proportion of links for each hub 
    #W = gamodels.RandomGraph(num_nrns, density*(num_nrns*(num_nrns-1)), directed=True)
    # just build the binary graph: the weights are added internally by the simulation program
    if is_numlike(seed):
        np.random.seed(seed)

    W = np.zeros([num_nrns, num_nrns])
    c = 0
    changeable_nodes = np.array(range(num_nrns))
    for i in hubs_list:
	#how many links you have to create according to hubs_density
        desired_links = np.round(hubs_density[c]*num_nrns)
	# how many links are already there
        actual_links = sum(W[:,i]==1)
        # random links to create 
        rnd_links = rnd.sample( changeable_nodes[W[changeable_nodes,i]==0], int(desired_links-actual_links) )
        # create links from incoming connections 
        W[rnd_links, i] = 1
	# update the list of nodes whose links can be changed
        changeable_nodes = changeable_nodes[changeable_nodes!=i]
	# how many links are already there
        actual_links = sum(W[i,:]==1)
        # random links to create
        rnd_links = rnd.sample( changeable_nodes[W[i,changeable_nodes]==0], int(desired_links-actual_links) )
        # create links from outgoing connections
        W[i, rnd_links] = 1
        c += 1
    # find changeable links
    mask = np.ones(np.shape(W),dtype=bool)
    for x in hubs_list:
	    mask[x,:] = False
	    mask[:,x] = False
    changeable_links = find(mask)  # linear indices of links non involved in hubs
    # flatten W
    tmp = np.reshape(W,-1)
    # how many links remain to create
    links_rem = size(W)*density - sum(tmp==1)
    # select random links to create
    rnd_links = rnd.sample( changeable_links, int(links_rem) )
    tmp[rnd_links] = 1
    W = reshape(tmp, shape(W))
    # reset the seed to a random value
    np.random.seed()
    return W

num_nrns = args.num_neurons 
density = args.density 
hubs = args.hubs_list
print(hubs)
hubs_density = args.hubs_density 
seed = 1

W = create_W(num_nrns, density, hubs, hubs_density, seed)

f_name = args.out_file 
f = gzip.open(f_name, 'wb')
pickle.dump(W, f, 1)
f.close()

[inD, outD] = galib.Degree(W, directed=True)
print('mean in degree:')
print(np.mean(inD))
print('mean out degree:')
print(np.mean(outD))
print('density:')
print(galib.Density(W))
figure
subplot(1,2,1)
hist(inD,60)
subplot(1,2,2)
hist(outD,60)
figure()
imshow(W, interpolation='nearest')
colorbar()
show()
