import numpy as np
import glob
import os
import argparse
import re
import networkx as nx
import scipy.io as sio
#import igraph as ig

def read_network(netfile):
    if '.mat' in netfile:
        # read network from matlab format file
        adjmat = sio.loadmat(netfile)['A']
        G = nx.from_numpy_matrix(adjmat)
    else:
        G = nx.read_edgelist(netfile, delimiter='\t', comments='#', data=(('weight', float),))
    return G

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname', '-d', type=str, required=True)
    parser.add_argument('--expname', '-e', type=str, default='exp_20180628')
    parser.add_argument('--infolder', '-i', type=str, default='mat')
    parser.add_argument('--outfolder', '-o', type=str, default='graphml')
    args = parser.parse_args()
    print(args)

    datname, expname, infolder, outfolder = args.dataname, args.expname, args.infolder, args.outfolder

    hfolder = '{}/{}'.format(datname, expname)
    infolder = '{}/{}'.format(hfolder, infolder)
    if not os.path.isdir(infolder):
        print('Input not found = {}'.format(infolder))
        exit(1)

    outfolder = '{}/{}'.format(hfolder, outfolder)
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    idx = 0
    for matfile in glob.glob(r"{}/nodes_*_index_*.*".format(infolder)):
        g = read_network(matfile)
        # convert via edge list
        #g1 = ig.Graph(len(g), zip(*zip(*nx.to_edgelist(g))[:2]))
        # nx.to_edgelist(g) returns [(0, 1, {}), (0, 2, {}), ...], which is turned
        #  into [(0, 1), (0, 2), ...] for igraph

        # convert via adjacency matrix
        #g2 = ig.Graph.Adjacency((nx.to_numpy_matrix(g) > 0).tolist())

        # write graph
        outfile = os.path.join(outfolder, os.path.basename(matfile))
        outfile = outfile.replace('.nse', '.graphml')
        outfile = outfile.replace('.mat', '.graphml')
        outfile = outfile.replace('_adj', '')
        outfile = outfile.replace('_avgdeg_16_maxdeg_32_', '_')
        if os.path.isfile(outfile):
            continue
        nx.write_graphml(g, outfile)
        idx = idx + 1
    print('Done with {} files'.format(idx))
