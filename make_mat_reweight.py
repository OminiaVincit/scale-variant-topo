import numpy as np
import glob
import os
import argparse
import re
import scipy.io as sio
from collections import defaultdict
from enum import Enum

from grakel import datasets


class REWEIGHT(Enum):
    NORMAL=0
    NEIGHBOR_LABELS_HIST=1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname', '-d', type=str, required=True)
    parser.add_argument('--exppath', '-e', type=str, default='exp_20190516')
    parser.add_argument('--outfolder', '-o', type=str, default='nse')
    parser.add_argument('--rw', type=str, default='nb_hist')
    parser.add_argument('--formattype', '-f', type=str, default='nse')
    args = parser.parse_args()
    print(args)

    exppath, datname, outfolder = args.exppath, args.dataname, args.outfolder
    ftype, rew = args.formattype, args.rw

    rw_lb = REWEIGHT.NORMAL
    if rew == 'nb_hist':
        rw_lb = REWEIGHT.NEIGHBOR_LABELS_HIST
    else:
        print('Reweight method is not defined. Keep weights as the normal')
    
    outfolder = '{}/{}/{}_{}'.format(exppath, datname, outfolder, rew)
    if not os.path.isdir(outfolder):
       os.makedirs(outfolder)

    dat = datasets.fetch_dataset(datname, as_graphs=True)
    G, y = dat.data, dat.target
    print(datname, y.shape)

    # find max node labels
    max_nlb = 0
    for i in range(y.shape[0]):
        lbgr = G[i].get_label_group(label_type="vertex", purpose="dictionary")
        max_nlb = max(max_nlb, max(lbgr.keys()))

    print('max_nlb = ', max_nlb)

    lbs = set(y)
    lbs = sorted(list(lbs))

    lbgr = G[0].get_label_group(label_type="vertex", purpose="dictionary")
    print(lbgr)
    print(G[0].get_labels(label_type="vertex", purpose="dictionary"))
    print(G[0].get_edges())
    print(G[0].get_vertices())


    lbgr = G[1].get_label_group(label_type="vertex", purpose="dictionary")
    print(lbgr)
    print(G[1].get_labels(label_type="vertex", purpose="dictionary"))
    print(G[1].get_edges())
    print(G[1].get_vertices())

    # get neighbor histogram for node labels
    for i in range(y.shape[0]):
        graph_id = i
        graph_lb = lbs.index(y[i])

        g = G[i]
        vlabels = g.get_labels(label_type="vertex", purpose="dictionary")
        voffset = min(vlabels.keys())

        edges = g.get_edges()
        vertices = g.get_vertices()
        
        # get list of labels of neighbor nodes for each node
        nbls = defaultdict(list)
        for v in vertices:
            nbls[v] = list()
        
        for e in edges:
            vfrom, vto = e[0], e[1]
            nbls[vfrom].append(vlabels[vto + voffset])    

        # get histogram of labels of neighbor nodes for each node
        nbhist = defaultdict(np.array)
        for v, nvl in nbls.items():
            tmp = np.array([0 for x in range(max_nlb+1)], dtype=np.float32)
            for ne_v in nvl:
               tmp[ne_v] += 1
            nbhist[v] = tmp
        # print(" ")
        # print(nbhist)
        # print(" ")

        # make distance between two adjacency nodes
        rs = []
        for e in edges:
            vfrom, vto = e[0], e[1]
            if vfrom < vto: # deal with undirected network
                w = 1.0
                if rw_lb == REWEIGHT.NEIGHBOR_LABELS_HIST:
                    diff = np.linalg.norm(nbhist[vfrom]-nbhist[vto])
                    w = np.exp(-np.square(diff)/2)
                rs.append('{}\t{}\t{}\n'.format(vfrom, vto, w))
        
        # write to file
        outfile = os.path.join(outfolder, 'nodes_{}_edges_{}_gid_{}_lb_{}_index_1_adj.nse'.format(len(vertices), len(edges)/2, graph_id, graph_lb))
        with open(outfile, 'w') as wf:
            wf.writelines(rs)

    
