######################################################################################################
## generate_diffusion_mat.py:  generate diffusion distance matrix from network adjency matrix
## Usage: python generate_diffusion_mat.py [options]
##      Options: --matpath : input directory consisting of the collection of adjency matrices 
#                   (nodes_*_index_{}_adj.mat --> for ER, GN, SP, WS, Drosophila-Melanogaster networks 
#                   nodes_*_index_{}_adj.nse  --> for LFR, LFR-H networks )
##              --outpath : output directory to save diffusion matrices
##                  (*.npy file with matrix size of K x N x N, K: number of tau, N: number of nodes)
##              --normflag : normalize (=1) the diffusion matrix by its maximum element or not
##              --avgflag  : average (=1) the diffusion matrix for K values of tau or not
##              --idxbg    : start index of adjency matrices file to divide in multi processing
##              --idxed    : end index of adjency matrices file to divide in multi processing
##              --taubg    : start value of tau
##              --taued    : end value of tau
##              --interval : sampling interval for tau from taubg to taued
##              --nproc    : number of processes using in this program 
##      Example: python generate_diffusion_mat.py 
#                           --matpath get-networks\mat\GN-net --outpath distmat\GN-net 
##                          --normflag 0  --avgflag 0  --idxbg 1      --idxed 1
##                          --taubg 1.0   --taued 10.0 --interval 1.0 --nproc 8
##
######################################################################################################

import networkx as nx
import numpy as np
import scipy.io as sio
import glob
import os
import argparse
from threading import Thread
from multiprocessing import Process

def read_network(netfile):
    if '.mat' in netfile:
        # read network from matlab format file
        adjmat = sio.loadmat(netfile)['A']
        G = nx.from_numpy_matrix(adjmat)
    else:
        G = nx.read_edgelist(netfile, delimiter='\t', comments='#', data=(('weight', float),))
    return G

def eig_decompose_mat(A):
    if A.shape[0] != A.shape[1]:
        raise Exception("Error: A is not square matrix!")
    
    if (A == A.T).all() == False:
        raise Exception("Error: A is not symmetric matrix!")
    
    la, U = np.linalg.eigh(A)
    U, _ = np.linalg.qr(U)
    if np.iscomplexobj(la):
        raise Exception("Error: la is not real matrix!")
    if np.iscomplexobj(U):
        raise Exception("Error after QR: U is not real matrix!")
    return U, la

def distmat(nzero_ls, zero_ls, U, la, Dleft, Dright, tau, normflag):
    D = np.diag(np.exp(-tau*la))
    Uleft = Dleft.dot(U)
    Uright = U.T.dot(Dright)
    Vsmall = Uleft.dot(D.dot(Uright))

    nz = len(nzero_ls) + len(zero_ls)
    V = np.identity(nz)
    V[np.ix_(nzero_ls, nzero_ls)] = Vsmall

    nr, nc = V.shape
    if nr != nc:
        raise Exception("Error: V is not square matrix")

    for i in range(nr):
        V[i] = np.abs(V[i])
        V[i] = V[i] / np.sum(V[i])
    
    Omega = np.zeros((nr, nc))
    for i in range(nr):
        for j in range(i, nr):
            Omega[i, j] = np.linalg.norm(V[i].flatten() - V[j].flatten())
    
    for i in range(nr):
        for j in range(0, i):
            Omega[i, j] = Omega[j, i]

    alpha = np.max(Omega)
    if normflag and alpha > 0:
        Omega = Omega / alpha
    return Omega

def distmat_taus(nzero_ls, zero_ls, Lsym, Dleft, Dright, tauls, normflag):
    U, la = eig_decompose_mat(Lsym)
    matls = []
    nz = len(nzero_ls) + len(zero_ls)
    for tau in tauls:
        Omega = distmat(nzero_ls, zero_ls, U, la, Dleft, Dright, tau, normflag)
        matls.append(Omega)
    matls = np.array(matls, dtype=np.float32).reshape(len(tauls), nz, nz)
    return matls

def tda_net(netfile, out_path, tauls, normflag, avgflag):
    basename = os.path.basename(netfile).replace('_adj.mat', '')
    basename = basename.replace('_adj.nse', '')
    basename = basename.replace('.nse', '')
    distpath = os.path.join(out_path, basename)

    G = read_network(netfile)
    A = nx.to_scipy_sparse_matrix(G)
    A = A.todense()

    nzero_ls = []
    zero_ls = []
    B = (A==0).all(1)
    nzero_ls = np.where(B==False)[0]
    zero_ls  = np.where(B==True)[0]
    
    if len(zero_ls) > 0:
        print('Need shrink the adjacency matrix then padding the probability matrix')
        print('netfile={}'.format(netfile))
        print('zero indexs', zero_ls)
    
    A = A[np.ix_(nzero_ls, nzero_ls)]
    diags = A.sum(axis=1).flatten()
    diagls = np.array(diags).ravel()
    a = np.sum(diagls==0)

    if a > 0:
        raise Exception("Error: There's an obsoleted node in network")
    
    D     = np.diag(diagls)
    Dleft = np.diag(1.0/np.sqrt(diagls))
    Dright = np.diag(np.sqrt(diagls))
    L = D - A
    Lsym = Dleft.dot(L.dot(Dleft))

    # make symmetric matrix for precision error
    Lsym = (Lsym + Lsym.T)/2.0

    matls = distmat_taus(nzero_ls, zero_ls, Lsym, Dleft, Dright, tauls, normflag)
    if avgflag:
        matls = np.mean(matls, axis=0, keepdims=True)
        print('Average of diffusion matrix, avg shape = ', matls.shape)
    np.save(distpath, matls)
    print('Saved {}'.format(os.path.basename(distpath)))

def get_file_list_process(mat_path, idxls):
    netfile_ls = []
    for idx in idxls:
        for netfile in glob.glob(r'{}/nodes_*_index_{}_adj.mat'.format(mat_path, idx)):
                # for ER, GN, SP, WS, Drosophila-Melanogaster networks 
                netfile_ls.append(netfile)
        for netfile in glob.glob(r'{}/nodes_*_index_{}_adj.nse'.format(mat_path, idx)):
                # for LFR, LFR-H networks
                netfile_ls.append(netfile)
    return netfile_ls

def make_exp(netfile_ls, out_path, tauls, normflag, avgflag, num_processes, process_id):
    print('Start proc id = {} per {} procs'.format(process_id, num_processes))
    if process_id >= num_processes:
        raise Exception('Proc id must be less than number of procs')
    for i in range(len(netfile_ls)):
        if i % num_processes == process_id:
            netfile = netfile_ls[i]
            print('Proc id={}'.format(process_id), i, os.path.basename(netfile))
            tda_net(netfile, out_path, tauls, normflag, avgflag)
    print('Finished proc id = {} per {} procs'.format(process_id, num_processes))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--matpath', type=str, required=True)
    parser.add_argument('--outpath', type=str, required=True)
    parser.add_argument('--normflag', type=int, default=0)
    parser.add_argument('--avgflag', type=int, default=0)
    parser.add_argument('--idxbg', type=int, default=1)
    parser.add_argument('--idxed', type=int, default=10)
    parser.add_argument('--taubg', type=float, required=True)
    parser.add_argument('--taued', type=float, required=True)
    parser.add_argument('--interval', type=float, required=True)
    parser.add_argument('--nproc', type=int, default=8)
    args = parser.parse_args()
    print(args)

    idx_bg, idx_ed = args.idxbg, args.idxed
    tau_bg, tau_ed, interval = args.taubg, args.taued, args.interval
    mat_path, out_path = args.matpath, args.outpath
    norm_flag, avg_flag, nproc = args.normflag, args.avgflag, args.nproc

    tauls = np.arange(tau_bg, tau_ed+interval, interval)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    np.savetxt(r'{}\norm_{}_avg_{}_ibg_{}_ied_{}_taubg_{}_taued_{}_interval_{}.txt'.format(out_path, norm_flag, avg_flag,
        idx_bg, idx_ed, tau_bg, tau_ed, interval), tauls)
    netfile_ls = get_file_list_process(mat_path, range(idx_bg, idx_ed+1))
    print('Number of file', len(netfile_ls))
    #
    # NOTICE: multi-thread is not scaled in python interpreter
    # See: https://medium.com/practo-engineering/threading-vs-multiprocessing-in-python-7b57f224eadb
    # for thread_id in range(nthreads):
        # p = threading.Thread(target=make_exp, args=(netfile_ls, out_path, tauls, nthreads, thread_id))
        # p.start()
    
    # LET'S USE Multi-processing
    processes = []
    for proc_id in range(0, nproc):
        p = Process(target=make_exp, args=(netfile_ls, out_path, tauls, norm_flag, avg_flag, nproc, proc_id))
        processes.append(p)
    
    # Start the processes
    for p in processes:
        p.start()
    
    # Ensure all processes have finished execution
    for p in processes:
        p.join()
        