import numpy as np
import glob
import os
import argparse
import re
import scipy.io as sio

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname', '-d', type=str, required=True)
    parser.add_argument('--expfolder', '-e', type=str, default='exp_20180628')
    parser.add_argument('--outfolder', '-o', type=str, default='mat')
    args = parser.parse_args()
    print(args)

    outfolder, expfolder, datname = args.outfolder, args.expfolder, args.dataname

    hfolder = '{}/{}'.format(datname, expfolder)
    outfolder = '{}/{}'.format(hfolder, outfolder)
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    adj_file = '{}/source/{}_A.txt'.format(hfolder, datname)
    indicator_file = '{}/source/{}_graph_indicator.txt'.format(hfolder, datname)
    graph_label_file = '{}/source/{}_graph_labels.txt'.format(hfolder, datname)

    # open graph label
    glabel = {}
    lbs = set()
    with open(graph_label_file, 'r') as grf:
        idx = 0
        for line in grf:
            idx = idx + 1
            line = line.strip()
            ar = [int(x) for x in line.split()]
            glabel[idx] = int(ar[0])
            lbs.add(glabel[idx])
    lbs = sorted(list(lbs))

    num_graphs = len(glabel.keys())
    gmap = {}
    netdict = {}

    for i in range(1, num_graphs+1):
        netdict[i] = []
        gmap[i] = []
    
    nodels = []
    # open graph indicator
    with open(indicator_file, 'r') as irf:
        node_id = 0
        for line in irf:
            node_id = node_id + 1
            line = line.strip()
            ar = [int(x) for x in line.split()]
            gid = ar[0]
            if gid > 0:
                gmap[gid].append(node_id)
                nodels.append(gid)
    
    # open adj file
    A = {}
    for gid in gmap.keys():
        N = len(gmap[gid])
        if N > 0:
            A[gid] = np.zeros((N, N), dtype=np.uint8)

    w = 1
    with open(adj_file, 'r') as arf:
        for line in arf:
            line = line.strip()
            line = line.replace(',', ' ')
            ar = [int(x) for x in line.split()]
            v1, v2 = ar[0], ar[1]
            gid = nodels[v1 - 1]
            n1 = v1 - min(gmap[gid])
            n2 = v2 - min(gmap[gid])
            A[gid][n1, n2] = 1
            A[gid][n2, n1] = 1
            netdict[gid].append('{}\t{}\t{}\n'.format(n1, n2, w))

    # save to mat file
    if True:
        for gid in A.keys():
            N = len(gmap[gid])
            ne = np.sum(A[gid])
            lb = glabel[gid]
            lb = lbs.index(lb)
            outfile = os.path.join(outfolder, 'nodes_{}_edges_{}_gid_{}_lb_{}_index_1_adj.mat'.format(N, ne, gid, lb))
            sio.savemat(outfile, {'A': A[gid]})
    
    # save to nse file
    else:
        for gid in netdict.keys():
            rs = netdict[gid]
            ne = len(rs)
            if ne > 0 :
                # save to file
                N = len(gmap[gid])
                lb = glabel[gid]
                outfile = os.path.join(outfolder, 'nodes_{}_edges_{}_gid_{}_lb_{}_index_1.nse'.format(N, ne, gid, lb))
                with open(outfile, 'w') as wf:
                    wf.writelines(rs)