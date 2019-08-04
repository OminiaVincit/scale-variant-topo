######################################################################################################
##
## make_common_features.py:  calculate common features from network adjency matrix
##      avg_deg_centrality, avg_clustering_coeff, avg_shortest_paths, 
##      avg_eccentricity, avg_node_closeness_centrality, 
##      avg_node_betweeness_centrality, avg_edge_betweeness_centrality, 
##      avg_local_efficiency, avg_triangles,
##      density, transitivity, diameter, radius, 
##      assortiative_coef, global_efficiency, nums_connected_parts
##
## Usage: python make_common_features.py -d dataname -i infolder -o outfolder
##      Example: python make_common_features.py -d MUTAG -i real-networks -o common-features
##
######################################################################################################

import numpy as np
import glob
import os
import argparse
import re
import scipy.io as sio
import community as cm #for modularity calculation

import pykov #for expectation of global mean first passage time
# https://github.com/riccardoscalco/Pykov

import networkx as nx
from multiprocessing import Process, Array
import tempfile
import time

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


def read_setting(filename):
    """
    Reading setting file by each type of network
    """
    namelist = []
    with open(filename, 'r', encoding='utf-8') as rf:
        lines = rf.readlines()
        for line in lines:
            line = line.strip()
            name = line.replace('graphml', 'mat')
            if '_adj.mat' not in name:
                name = name.replace('.mat', '_adj.mat')
            namelist.append(name)
        rf.close()
    print('Number of files = {}'.format(len(namelist)))
    return namelist

def read_network(netfile):
    if '.mat' in netfile:
        # read network from matlab format file
        adjmat = sio.loadmat(netfile)['A']
        G = nx.from_numpy_matrix(adjmat)
    else:
        G = nx.read_edgelist(netfile, delimiter='\t', comments='#', data=(('weight', float),))
    return G

def sum_dict(dic):
    total = 0.0
    for a in dic:
        total += dic[a]
    return total

def get_modularity(G):
    """
    Compute the maximum modularity using the Louvain heuristices
    obtain different results at each run
    """
    part = cm.best_partition(G)
    md = cm.modularity(part, G)
    return md

def get_EGMFPT(G):
    """
    Get expectation of global mean first passage time
    """
    A = nx.to_scipy_sparse_matrix(G)
    A = A.todense()
    H, L = A.shape[0], A.shape[1]
    #P = np.zeros((H, L), dtype=np.float32)
    dtrans = {}
    for i in range(H):
        ki = np.sum(A[i, :])
        for j in range(L):
            if ki == 0 or i == j:
                dtrans[(i, j)] = 0
            else:
                dtrans[(i, j)] = A[i, j] / ki
    T = pykov.Chain(dtrans)
    mfpts = [np.mean([v for a, v in T.mfpt_to(i).items()]) for i in range(H)]
    egm = np.average(mfpts)
    if np.isnan(egm) or egm < 0 or egm >= 1e+6:
        print('Mean of global mean first passage time not defined')
        egm = 0.0
    # for i in range(H):
    #     print(i, T.mfpt_to(i))
    return egm

def network_statistics(G):
    # Average of local features
    num_nodes = len(G)
    netsize = G.size()
    features = []
    
    avg_deg_centrality = float(sum_dict(nx.degree_centrality(G))) / float(num_nodes)
    avg_cc = nx.average_clustering(G)
    try:
        avg_shortest = nx.average_shortest_path_length(G)
    except nx.NetworkXError:
        avg_shortest = 0
    
    try:
        avg_eccentricity = float(sum_dict(nx.eccentricity(G)))/float(num_nodes)
    except nx.NetworkXError:
        avg_eccentricity = 0.0
    
    avg_node_closeness_centrality = float(sum_dict(nx.closeness_centrality(G)))/float(num_nodes)
    avg_node_betweeness_centrality = float(sum_dict(nx.betweenness_centrality(G)))/float(num_nodes)
    avg_edge_betweeness_centrality = float(sum_dict(nx.edge_betweenness_centrality(G)))/float(netsize)
    avg_local_efficiency = nx.local_efficiency(G)
    avg_triangles = float(sum_dict(nx.triangles(G)))/float(num_nodes)
    
    # Global features
    density  = float(2*netsize)/float(num_nodes*(num_nodes-1))
    transitivity = nx.transitivity(G)
    try:
        diameter = nx.diameter(G)
        radius = nx.radius(G)
    except nx.NetworkXError:
        diameter = 0
        radius = 0
    
    ass_coef = nx.degree_assortativity_coefficient(G)
    global_efficiency = nx.global_efficiency(G)
    nums_connected = nx.number_connected_components(G)
    md = get_modularity(G)
    #egm = get_EGMFPT(G)
    egm = 0.0
    features = [egm, md, ass_coef, avg_cc, avg_deg_centrality, avg_shortest, avg_eccentricity, avg_node_closeness_centrality, 
                avg_node_betweeness_centrality, avg_edge_betweeness_centrality, avg_local_efficiency, avg_triangles,
                density, transitivity, diameter, radius, global_efficiency, nums_connected]
    return features 


def save_features_process(filels, a, shape, bg_idx, process_id):
    print('Enter process {}'.format(process_id))
    a = np.frombuffer(a.get_obj(), dtype=np.float32)
    a.shape = shape
    features = []
    count = 1
    for netfile in filels:
        G = read_network(netfile)
        features.append(network_statistics(G))
        if count % 10 == 0:
            print('Process {} has finished to process {} / {} files'.format(process_id, count, len(filels)))
        count += 1
    a[bg_idx:(bg_idx + len(filels))] = np.asarray(features, dtype=np.float32)
    print('Finished process {} begin at {} with {} files'.format(process_id, bg_idx, len(filels)))
    # tmpfile = '{}/features_{}'.format(tmpfolder, process_id)
    # np.save(tmpfile, arr)
    # print('Saved temporarily {}'.format(tmpfile))


def make_features(datname, infolder, outfolder, n_jobs=-1):
    start = time.time()
    graphlist = os.path.join(infolder, '{}/graphlist.txt'.format(datname))

    if not os.path.isfile(graphlist):
        # if not os.path.isfile(graphlist):
        #     graphlist = os.path.join(infolder, '{}/matlist_grakel.txt'.format(datname))
        if not os.path.isfile(graphlist):
            print('File not found: {}'.format(graphlist))

    namels = read_setting(graphlist)
    N = len(namels)
    features = []
    if N > 0 and os.path.isfile(namels[0]):
        filels = namels
    else:
        filels = ['{}/{}/{}'.format(infolder, datname, name) for name in namels]

    # 18 features for each network
    arr = np.zeros((N, 18), dtype=np.float32)

    if n_jobs <= 0:
        for i in range(N):
            netfile = filels[i]
            G = read_network(netfile)
            features.append(network_statistics(G))
            print('Processed {} files = {}'.format(i+1, netfile))
        arr[...] = np.asarray(features)

    else:
        # Multiprocessing version
        split_ls = [list(b) for b in np.array_split(range(N), n_jobs)]
        split_ls = [[filels[i] for i in b] for b in split_ls]
        shape = arr.shape
        size = arr.size
        arr.shape = size
        marr = Array('f', arr)
        #with tempfile.TemporaryDirectory() as temp_path:

        processes = []
        bg_idx = 0
        for proc_id in range(0, n_jobs):
            p = Process(target=save_features_process, args=(split_ls[proc_id], marr, shape, bg_idx, proc_id))
            processes.append(p)
            bg_idx += len(split_ls[proc_id])
        
        # Start the processes
        for p in processes:
            p.start()
        
        # Ensure all processes have finished execution
        for p in processes:
            p.join()

        marr = np.frombuffer(marr.get_obj(), dtype=np.float32)
        marr.shape = shape
        arr.shape = shape
        arr[...] = marr
        
        # Load all tempfile
        # for proc_id in range(0, n_jobs):
        #     fet_arr = np.load('{}/features_{}.npy'.format(temp_path, proc_id), )
        #     features.extend([list(a) for a in fet_arr])

    print('Feature mat', arr.shape)

    outpath = '{}'.format(outfolder)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    features_file = '{}/{}_common_features_rev'.format(outpath, datname)
    np.save(features_file, arr)
    elapse = sec_to_time(round(time.time() - start, 2))
    print('Calculated in {} to save {} '.format(elapse, features_file))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname',  '-d', type=str, required=True)
    parser.add_argument('--outfolder', '-o', type=str, default='common-features')
    parser.add_argument('--infolder',  '-i', type=str, default='real-networks')
    parser.add_argument('--njobs',  '-n', type=int, default=0)
    args = parser.parse_args()
    print(args)

    outfolder, infolder, datname = args.outfolder, args.infolder, args.dataname
    #G = nx.Graph()
    #G.add_edges_from([(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (3, 4), (4, 3), (5, 6)])
    #G.add_node(5)
    #get_EGMFPT(G)
    make_features(datname, infolder, outfolder, args.njobs)