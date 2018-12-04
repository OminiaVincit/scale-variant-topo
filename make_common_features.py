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
import networkx as nx

def read_setting(filename):
    """
    Reading setting file by each type of network
    """
    namelist = []
    with open(filename, 'r') as rf:
        lines = rf.readlines()
        for line in lines:
            line = line.strip()
            name = line.replace('graphml', 'mat')
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
    
    features = [avg_deg_centrality, avg_cc, avg_shortest, avg_eccentricity, avg_node_closeness_centrality, 
                avg_node_betweeness_centrality, avg_edge_betweeness_centrality, avg_local_efficiency, avg_triangles,
                density, transitivity, diameter, radius, ass_coef, global_efficiency, nums_connected]
    return features 

def make_features(datname, infolder, outfolder):
    graphlist = os.path.join(infolder, '{}/graphlist.txt'.format(datname))

    if not os.path.isfile(graphlist):
        print('File not found: {}'.format(graphlist))

    namels = read_setting(graphlist)
    N = len(namels)
    features = []
    
    for i in range(N):
        netfile = '{}/{}/{}'.format(infolder, datname, namels[i])
        G = read_network(netfile)
        features.append(network_statistics(G))
        print('Processed {} files = {}'.format(i+1, netfile))

    arr = np.asarray(features)
    print('Feature mat', arr.shape)

    outpath = '{}'.format(outfolder)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    features_file = '{}/{}_common_features'.format(outpath, datname)
    np.save(features_file, arr)
    print('Saved {}'.format(features_file))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname',  '-d', type=str, required=True)
    parser.add_argument('--outfolder', '-o', type=str, default='common-features')
    parser.add_argument('--infolder',  '-i', type=str, default='real-networks')
    args = parser.parse_args()
    print(args)

    outfolder, infolder, datname = args.outfolder, args.infolder, args.dataname
    make_features(datname, infolder, outfolder)