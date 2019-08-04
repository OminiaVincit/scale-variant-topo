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

from collections import defaultdict
import argparse

#from multiprocessing import shared_memory

DS_LB, CON_LB = 'discrete', 'continous'

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

def node_features_from_graph(G, labels, prefer_attr_nodes):
    matdict = G.get_labels()
    gls = [ v for s, v in sorted(matdict.items()) ]
    if prefer_attr_nodes == False:
        n_labels = len(labels)
        #print('nlabels', n_labels, gls)
        gls = [ np.eye(n_labels)[labels.index(v)] for v in gls ]
    n = len(gls)
    if n > 0:
        m = len(gls[0])
        gls = np.asarray(gls, dtype=np.float32).reshape(n, m)
    return gls

def graph_features_from_dataset(dataname, prefer_attr_nodes):
    gdat = datasets.fetch_dataset(dataname, prefer_attr_nodes=prefer_attr_nodes, as_graphs=True).data
    labels = []
    if prefer_attr_nodes == False:
        lbs = set()
        for G in gdat:
            for lb in G.get_label_group().keys():
                lbs.add(lb)
        labels = sorted(list(lbs))
    rs = defaultdict(int)
    for i in range(len(gdat)):
        rs[i] = node_features_from_graph(gdat[i], labels=labels, prefer_attr_nodes=prefer_attr_nodes)
    return rs

def discrete_kernel(mat_a, mat_b):
    return np.matmul(mat_a, mat_b.T)

def pairwise_continous_kernel(vec_a, vec_b, p, gamma):
    return np.exp(-gamma*np.linalg.norm(vec_a-vec_b, ord=p))

def continous_kernel(mat_a, mat_b, p, gamma):
    nsz_a, nsz_b = mat_a.shape[0], mat_b.shape[0]
    K = np.zeros((nsz_a, nsz_b), dtype=np.float32)
    for i in range(nsz_a):
        for j in range(nsz_b):
            K[i, j] = pairwise_continous_kernel(mat_a[i], mat_b[j], p, gamma)
    return K


def pairwise_domains_graph_kernels(mat_a, mat_b, kertype, p, gamma):
    K = None
    if kertype == 0:
        # only discrete kernel for label
        K = discrete_kernel(mat_a[DS_LB], mat_b[DS_LB])
    elif kertype == 1:
        # only continous kernel for attributes
        K = continous_kernel(mat_a[CON_LB], mat_b[CON_LB], p, gamma)
    else:
        # tensor of discrete kernel and continous kernel
        K_d = discrete_kernel(mat_a[DS_LB], mat_b[DS_LB])
        K_c = continous_kernel(mat_a[CON_LB], mat_b[CON_LB], p, gamma)
        K = np.multiply(K_d, K_c)
    return K

def calculate_graph_embedd_distance(glist, kertype, p, gamma):
    start = time.time()
    N = len(glist)
    D = np.zeros((N, N), dtype=np.float32)

    for i in range(N):
        G = glist[i]
        K_gg = pairwise_domains_graph_kernels(G, G, kertype, p, gamma)
        N_g = K_gg.shape[0]
        one_g = np.ones(N_g)
        for j in range(i+1, N):
            H = glist[j]
            K_hh = pairwise_domains_graph_kernels(H, H, kertype, p, gamma)
            K_gh = pairwise_domains_graph_kernels(G, H, kertype, p, gamma)
            
            N_h = K_hh.shape[0]
            one_h = np.ones(N_h)

            D[i, j] = np.sqrt(np.dot(np.dot(one_g, K_gg), one_g) / (N_g * N_g) + np.dot(np.dot(one_h, K_hh), one_h) / (N_h * N_h) \
                - 2 * np.dot(np.dot(one_g, K_gh), one_h) / (N_g * N_h))
            D[j, i] = D[i, j]
        elapse = sec_to_time(round(time.time() - start, 2))
        print('Graph kernel i={}, j={}, N={}, {}'.format(i, j, N, elapse))
    return D

# def element_process(G, H, i, j, p, gamma, shared_D):


# def multi_process_graph_embedd_distance(glist, kertype, p, gamma):
#     start = time.time()
#     N = len(glist)
#     D = np.zeros((N, N), dtype=np.float32)

#     shm = shared_memory.SharedMemory(create=True, size = D.nbytes)
#     B = np.ndarray((N, N), dtype=D.dtype, buffer=shm.buf)
#     B[:] = D[:]


#datalist = ['MUTAG', 'BZR', 'COX2', 'DHFR', 'ENZYMES', 'PROTEINS', 'NCI1', 'NCI109', 'DD', 'MSRC_9']

def save_kernel(outfile, glist, kertype, p_local, gamma_local, p_global, gamma_global):
    D = calculate_graph_embedd_distance(glist, kertype, p_local, gamma_local)
    print('Begin saving kernel = {}'.format(outfile))
    Dp = np.power(D, p_global)
    
    if gamma_global < 0:
        med = np.median(Dp.flatten())
        if med > 0:
            gamma_global = 1.0 / med
            print('Set global gamma={} with med = {}'.format(gamma_global, med))
        else:
            gamma_global = 1.0
            print('Set default global gamma=1.0')
    K = np.exp(- gamma_global * Dp)
    print(K)
    np.savetxt(outfile, K)
    print('Saved kernel shape=', K.shape)
    return K

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exppath', '-e', type=str, required=True)
    parser.add_argument('--folder', '-f', type=str, default='gkernel')
    parser.add_argument('--njobs', '-j', type=int, default=-1)
    parser.add_argument('--p1', type=int, default=2)
    parser.add_argument('--g1', type=float, default=-1.0)
    parser.add_argument('--p2', type=int, default=2)
    parser.add_argument('--g2', type=float, default=-1.0)
    parser.add_argument('--kertype', '-k', type=int, default=0) #0: only discrete for labels, #1: only continous for attr, #2: both 
    parser.add_argument('--dataname', '-d', type=str, default='')
    
    args = parser.parse_args()
    print(args)
    dataname = args.dataname
    p1, g1, p2, g2, kertype = args.p1, args.g1, args.p2, args.g2, args.kertype

    outpath = os.path.join(args.exppath, dataname)
    outpath = os.path.join(outpath, args.folder)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    glist = defaultdict(int)

    kerstr = ''
    glbs, gats = None, None

    if kertype == 0:
        kerstr = 'label'
        glbs = graph_features_from_dataset(dataname, prefer_attr_nodes=False)
        for gid, fet in sorted(glbs.items()):
            glist[gid] = defaultdict()
            glist[gid][DS_LB] = fet

    elif kertype == 1:
        kerstr = 'attr'
        gats = graph_features_from_dataset(dataname, prefer_attr_nodes=True)
        for gid, fet in sorted(gats.items()):
            glist[gid] = defaultdict()
            glist[gid][CON_LB] = fet
    else:
        kerstr = 'label_attr'
        glbs = graph_features_from_dataset(dataname, prefer_attr_nodes=False)
        gats = graph_features_from_dataset(dataname, prefer_attr_nodes=True)

        for gid, fet in sorted(glbs.items()):
            glist[gid] = defaultdict()
            glist[gid][DS_LB] = fet
            glist[gid][CON_LB] = gats[gid]
    
    if g1 < 0 and gats is not None:
        g1 = np.sqrt(gats[0].shape[1])
        print('Set local gamma={}'.format(g1))
    
    outfile = os.path.join(outpath, '{}_Hilbert_{}_p1_{}_p2_{}.txt'.format(dataname, kerstr, p1, p2))
    save_kernel(outfile, glist, kertype, p1, g1, p2, g2)