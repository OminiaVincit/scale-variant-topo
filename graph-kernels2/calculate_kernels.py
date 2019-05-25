import numpy as np
import os

from time import time

from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score

from grakel import datasets
from grakel import GraphKernel
from grakel.kernels import VertexHistogram, ShortestPath, WeisfeilerLehman, RandomWalkLabeled

import argparse

kernels = {
    "GraphletSampling": [{"name": "graphlet_sampling", "sampling": {"n_samples": 150}}],
    "WL-Subtree": [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}],
    "WL-ShortestPath": [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "shortest_path"}]
}

#gk = WeisfeilerLehman(niter=1, normalize=True, base_kernel=VertexHistogram)
#gk = VertexHistogram(normalize=True)

def save_kernel(G, gk, outpath, dataname, kername, handle=False):
    print('Compute kernel ', kername)
    n = len(G)
    b = 188 
    # Let's use multi-processing
    processes = []

    if handle == True:
        K = np.zeros((n,n))
        for i in range(0, n, b):
            ib = min(n, i+b)
            Gs = G[i:ib]
            Ks = gk.fit_transform(Gs)
            K[i:ib, i:ib] = Ks
            Gn = G[ib:n]
            if len(Gn) > 0:
                Kn = gk.transform(Gn)
                K[i:ib, ib:n] = Kn.T
                K[ib:n, i:ib] = Kn
    else:
        K = gk.fit_transform(G)
    P = gk.fit_transform(G)
    print('K')
    print(K)
    print('P')
    print(P)
    print('K-P')
    print(np.abs(P-K))

    outfile = os.path.join(outpath, '{}_{}.txt'.format(dataname, kername))
    #np.savetxt(outfile, K)
    #print('Saved kernel ', kername, K.shape)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exppath', '-e', type=str, required=True)
    parser.add_argument('--folder', '-f', type=str, default='gkernel')
    parser.add_argument('--njobs', '-j', type=int, default=8)
    args = parser.parse_args()
    print(args)
    njobs = args.njobs

    datalist = ['MUTAG', 'BZR', 'COX2', 'DHFR', 'ENZYMES', 'PROTEINS', 'NCI1', 'NCI109', 'DD', 'MSRC_9']
    datalist = ['NCI1', 'NCI109']
    datalist = ['DD']
    datalist = ['MUTAG']
    rows = sorted(list(kernels.keys()))

    for dataname in datalist:
        outpath = os.path.join(args.exppath,dataname)
        outpath = os.path.join(outpath, args.folder)
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
            
        dat = datasets.fetch_dataset(dataname, as_graphs=True)
        G, y = dat.data, dat.target
        print(dataname, y.shape)
        
        # Need to run each of below kernels separately
        if False:

            gk = VertexHistogram(normalize=True)
            save_kernel(G, gk, outpath, dataname, 'VertexHist', n_jobs=None)

            gk = ShortestPath(normalize=True)
            save_kernel(G, gk, outpath, dataname, 'ShortestPath', n_jobs=None)

        if False:
            gk = WeisfeilerLehman(niter=5, normalize=True, base_kernel=VertexHistogram, n_jobs=None)
            save_kernel(G, gk, outpath, dataname, 'WL-VertexHist')

        # if False:
        #     for rwtype in ['geometric', 'exponential']:
        #         gk = RandomWalkLabeled(normalize=True, kernel_type=rwtype)
        #         save_kernel(G, gk, outpath, dataname, 'randomwalk_{}'.format(rwtype))

        if True:
            for (i, kname) in enumerate(rows):
                print(kname, end=" ")
                gk = GraphKernel(kernel=kernels[kname], normalize=False, n_jobs=None)
                print("", end=".")
                save_kernel(G, gk, outpath, dataname, kname.replace('/', '-'), True)
                break
