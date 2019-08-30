######
### This script only works with the version of networkx < 2.0
######
import numpy as np
import glob
import os
import argparse
import re
import scipy.io as sio
import networkx as nx

from networkx.utils import powerlaw_sequence

spath = r'F:\Research\ScaleVariant\exp_20190803\scalefree'
N = 200
#model ='conf' 
model ='BA'

if model == 'BA':
  for a in range(1, 6):
    m = a
    for i in range(1, 21):
      filename = os.path.join(spath, r'{0}/nodes_{1}_{0}_m_{2}_index_{3}.txt'.format(model, N, m, i))
      G = nx.barabasi_albert_graph(N, m)
      nx.write_edgelist(G, filename, data=False)

if model == 'conf':
  for m in range(1, 101):        
    z = nx.utils.create_degree_sequence(N, powerlaw_sequence, exponent=3.0)
    G = nx.configuration_model(z)
    filename = os.path.join(spath, r'{0}/nodes_{1}_{0}_index_{2}.txt'.format(model, N, m))
    
    # remove parallel edges and self-loops
    loops = G.selfloop_edges()
    G = nx.Graph(G)
    G.remove_edges_from(loops)

    nx.write_edgelist(G, filename, data=False)