import numpy as np
import os

import time

from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score

from grakel import datasets
from grakel import GraphKernel
from grakel.kernels import VertexHistogram, ShortestPath, WeisfeilerLehman, RandomWalkLabeled, MultiscaleLaplacianFast

from six import itervalues, iteritems

import argparse

def sec_to_time(sec):
    """Print time in a correct format."""
    dt = list()
    days = int(sec // 86400)
    if days > 0:
        sec -= 86400*days
        dt.append(str(days) + " d")

    hrs = int(sec // 3600)
    if hrs > 0:
        sec -= 3600*hrs
        dt.append(str(hrs) + " h")

    mins = int(sec // 60)
    if mins > 0:
        sec -= 60*mins
        dt.append(str(mins) + " m")

    if sec > 0:
        dt.append(str(round(sec, 2)) + " s")
    return " ".join(dt)

lb_kernels = {
    "GraphletSampling": [{"name": "graphlet_sampling", "sampling": {"n_samples": 150}}],
    "WL-Subtree": [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}],
    "WL-ShortestPath": [{"name": "weisfeiler_lehman", "niter": 5}, {"name": "shortest_path"}]
}

ulb_kernels = {
    "ShortestPath" : [{"name": "shortest_path", "with_labels": False}],
    "GraphletSampling": [{"name": "graphlet_sampling", "sampling": {"n_samples": 150}}],
    "GeometricRandomWalk" : [{"name": "random_walk", "method_type": "fast", "with_labels": False, "kernel_type": "geometric"}], #ill defined, donot normalize
    "ExponentialRandomWalk" : [{"name": "random_walk", "method_type": "fast", "with_labels": False, "kernel_type": "exponential"}],
    # Must have node attribute "MultiScaleLaplacianFast" : [{"name": "multiscale_laplacian", "which": "fast"}],
    "LovaszTheta" : [{"name": "lovasz_theta"}], #slow
    #"SvmTheta" : [{"name": "svm_theta"}] #fast
}
#gk = WeisfeilerLehman(niter=1, normalize=True, base_kernel=VertexHistogram)
#gk = VertexHistogram(normalize=True)

def save_kernel(G, gk, outpath, dataname, kername, b, handle=False):
    start = time.time()
    print('Compute kernel {} use handle = {}'.format(kername, handle))
    n = len(G)
    #TODO: Let's use multi-processing but need to handle with large memory consuming problem
    K = np.zeros((n,n))
    if handle == True:
        for i in range(0, n, b):
            ib = min(n, i+b)
            Gs = G[i:ib]
            Ks = gk.fit_transform(Gs)
            K[i:ib, i:ib] = Ks
            for j in range(ib, n, b):
                jb = min(n, j+b)
                Gn = G[j:jb]
                Kn = gk.transform(Gn)
                K[i:ib, j:jb] = Kn.T
                K[j:jb, i:ib] = Kn
                elapse = sec_to_time(round(time.time()-start, 2))
                print('i={}, j={}, b={}, {}'.format(i, j, b, elapse))
    else:
        K = gk.fit_transform(G)
    # P = gk.fit_transform(G)
    # print('K')
    # print(K)
    # print('P')
    # print(P)
    # print('K-P')
    # print(np.max(np.abs(P-K)))

    outfile = os.path.join(outpath, '{}_{}.txt'.format(dataname, kername))
    end = time.time()
    elapse = sec_to_time(round(end-start, 2))
    print('Calculate kernel {} in {} '.format(kername, elapse))
    np.savetxt(outfile, K)
    print('Saved kernel ', kername, K.shape)
    print('')

def to_one_hot(G):
    # Index all discrete labels
    mp = {dl: i for (i, dl) in enumerate(set(l for g in G for l in itervalues(g[1])))}
    def make_vec(k):
        vec = np.zeros((len(mp),), dtype=float)
        vec[k] = 1.0
        return vec
    return [(g[0], {i: make_vec(mp[k]) for (i, k) in iteritems(g[1])}) for g in G]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exppath', '-e', type=str, required=True)
    parser.add_argument('--folder', '-f', type=str, default='gkernel')
    parser.add_argument('--njobs', '-j', type=int, default=-1)
    parser.add_argument('--norm', '-n', type=int, default=1)
    parser.add_argument('--handle', type=int, default=0)
    parser.add_argument('--batchsize', '-b', type=int, default=128)
    parser.add_argument('--label', type=int, default=0)
    parser.add_argument('--dataname', '-d', type=str, default='')
    parser.add_argument('--kername', '-k', type=str, default='')
    
    args = parser.parse_args()
    print(args)
    njobs = None
    norm, handle, b, label = args.norm, args.handle, args.batchsize, args.label
    if args.njobs > 0:
        njobs = args.njobs

    dname, kname = args.dataname, args.kername

    lb_datalist = ['MUTAG', 'BZR', 'COX2', 'DHFR', 'ENZYMES', 'PROTEINS', 'NCI1', 'NCI109', 'DD', 'MSRC_9']
    ulb_datalist = ['IMDB-BINARY', 'IMDB-MULTI', 'REDDIT-BINARY','FRANKENSTEIN', 'COLLAB']

    if label > 0:
        datalist = lb_datalist
        kernels = lb_kernels
    else:
        datalist = ulb_datalist
        kernels = ulb_kernels

    rows = sorted(list(kernels.keys()))

    if dname != '' and dname not in datalist:
        raise ValueError('Not found specified data: {}'.format(dname))
    
    if kname != '' and kname not in kernels:
        raise ValueError('Not found specified kernel: {}'.format(kname))

    for dataname in datalist:
        if dname != '' and dataname != dname:
            continue
        outpath = os.path.join(args.exppath,dataname)
        outpath = os.path.join(outpath, args.folder)
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
            
        dat = datasets.fetch_dataset(dataname, as_graphs=True)
        G, y = dat.data, dat.target
        print(dataname, y.shape)
        
        # Need to run each of below kernels separately
        if False and label > 0:
            gk = VertexHistogram(normalize=norm, n_jobs=njobs)
            save_kernel(G, gk, outpath, dataname, 'VertexHist', b, handle=handle)

            gk = ShortestPath(normalize=norm, n_jobs=njobs)
            save_kernel(G, gk, outpath, dataname, 'ShortestPath', b, handle=handle)

        if False and label > 0:
            gk = WeisfeilerLehman(niter=5, normalize=norm, base_kernel=VertexHistogram, n_jobs=None)
            save_kernel(G, gk, outpath, dataname, 'WL-VertexHist', b, handle=handle)

        # if False:
        #     for rwtype in ['geometric', 'exponential']:
        #         gk = RandomWalkLabeled(normalize=True, kernel_type=rwtype)
        #         save_kernel(G, gk, outpath, dataname, 'randomwalk_{}'.format(rwtype))

        if True:
            for (i, kername) in enumerate(rows):
                if kname != '' and kername != kname:
                    continue
                print(kername, end=" ")
                gk = GraphKernel(kernel=kernels[kername], normalize=norm, n_jobs=njobs)
                print("", end=".")
                use_handy = False
                if 'WL' in kername and len(G) > 256:
                    use_handy = True
                save_kernel(G, gk, outpath, dataname, kername.replace('/', '-'), b, handle=use_handy)
                
