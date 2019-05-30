import numpy as np
import glob
import os
import argparse
import re
from collections import defaultdict
from grakel import datasets
import scipy.io as sio

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exppath', '-e', type=str, required=True)
    parser.add_argument('--folder', '-f', type=str, default='grakel')
    parser.add_argument('--data', '-d', type=str, default='')
    args = parser.parse_args()
    print(args)

    labels_datals = ['MUTAG', 'BZR', 'COX2', 'DHFR', 'ENZYMES', 'DD', 'PROTEINS', 'NCI1', 'NCI109', 'MSRC_9']
    unlabels_datals = ['FRANKENSTEIN', 'IMDB-BINARY', 'IMDB-MULTI', 'REDDIT-BINARY', 'REDDIT-MULTI-5K', 'REDDIT-MULTI-12K']

    datalist = unlabels_datals
    
    for dataname in datalist:
        outpath = os.path.join(args.exppath, '{}/{}'.format(dataname, args.folder))

        if not os.path.isdir(outpath):
            os.mkdir(outpath)
            
        dat = datasets.fetch_dataset(dataname, as_graphs=True)
        G, y = dat.data, dat.target
        print(dataname, y.shape)

        lbs = set(y)
        lbs = sorted(list(lbs))

        for i in range(y.shape[0]):
            g = G[i]
            lb = lbs.index(y[i])
            A = g.get_adjacency_matrix()
            A = np.array(A, dtype=np.uint8)
            N = A.shape[0]
            ne = np.sum(A)

            outfile = os.path.join(outpath, 'nodes_{}_edges_{}_gid_{}_lb_{}_index_1_adj.mat'.format(N, ne, i, lb))
            sio.savemat(outfile, {'A': A})

