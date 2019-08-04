import numpy as np
import glob
import os
import argparse
import re
import scipy.io as sio
from collections import defaultdict
from enum import Enum

class REWEIGHT(Enum):
    NORMAL=0
    NEIGHBOR_LABELS_HIST=1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname', '-d', type=str, required=True)
    parser.add_argument('--exppath', '-e', type=str, default='exp_20180628')
    parser.add_argument('--outfolder', '-o', type=str, default='mat')
    parser.add_argument('--rw', type=str, default='normal')
    parser.add_argument('--formattype', '-f', type=str, default='mat')
    args = parser.parse_args()
    print(args)

    exppath, datname, outfolder = args.exppath, args.dataname, args.outfolder
    ftype, rew = args.formattype, args.rw

    rw_lb = REWEIGHT.NORMAL
    if rew == 'nb_hist':
        rw_lb = REWEIGHT.NEIGHBOR_LABELS_HIST
    else:
        print('Reweight method is not defined. Keep weights as the normal')
    
    outfolder = '{}/{}_{}'.format(exppath, outfolder, rew)
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    adj_file = '{}/source/{}_A.txt'.format(exppath, datname)
    indicator_file = '{}/source/{}_graph_indicator.txt'.format(exppath, datname)
    graph_label_file = '{}/source/{}_graph_labels.txt'.format(exppath, datname)
    node_label_file = '{}/source/{}_node_labels.txt'.format(exppath, datname)
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
    
    nodes_to_gid = []
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
                nodes_to_gid.append(gid)

    # open node label file
    nodes_to_label = []
    with open(node_label_file, 'r') as nrf:
        for line in nrf:
            line = line.strip()
            ar = [int(x) for x in line.split()]
            if len(ar) > 0:
                nodes_to_label.append(ar[0])
    max_nlb = max(nodes_to_label)

    # open adj file
    A = {}
    for gid in gmap.keys():
        N = len(gmap[gid])
        if N > 0:
            A[gid] = np.zeros((N, N), dtype=np.uint8)


    with open(adj_file, 'r') as arf:
        for line in arf:
            line = line.strip()
            line = line.replace(',', ' ')
            ar = [int(x) for x in line.split()]
            v1, v2 = ar[0], ar[1]
            gid = nodes_to_gid[v1 - 1]
            n1 = v1 - min(gmap[gid])
            n2 = v2 - min(gmap[gid])
            A[gid][n1, n2] = 1
            A[gid][n2, n1] = 1
            #netdict[gid].append('{}\t{}\t{}\n'.format(n1, n2, w))

    # save to mat file
    if ftype == 'mat':
        for gid in A.keys():
            N = len(gmap[gid])
            ne = np.sum(A[gid])
            lb = glabel[gid]
            lb = lbs.index(lb)
            outfile = os.path.join(outfolder, 'nodes_{}_edges_{}_gid_{}_lb_{}_index_1_adj.mat'.format(N, ne, gid, lb))
            sio.savemat(outfile, {'A': A[gid]})
    
    # save to nse file
    elif ftype == 'nse':
        count = 0
        for gid in A.keys():
            gmin = min(gmap[gid])
            B = np.transpose(np.nonzero(A[gid]))
            #print('Graph id ', gid, B)
            nb = defaultdict(np.array)
            deg = defaultdict(int)

            for i in range(A[gid].shape[0]):
                nb[i] = np.array([0 for x in range(max_nlb+1)], dtype=np.float32)
                deg[i] = 0
                #nb[i].append(nodes_to_label[i + gmin - 1])

            for b in B:
                i, j = b[0], b[1]
                if i < j:
                    lb_i, lb_j = nodes_to_label[i + gmin - 1], nodes_to_label[j + gmin - 1]
                    nb[i][lb_j] += 1
                    nb[j][lb_i] += 1
                    deg[i] += 1
                    deg[j] += 1
            
            # write to file
            rs = []
            for b in B:
                i, j = b[0], b[1]
                if i < j :
                    #nb[i] = nb[i]/np.sum(nb[i])
                    #nb[j] = nb[j]/np.sum(nb[j])
                    
                    #print('Node ', i+1, nb[i])
                    #print('Node ', j+1, nb[j])
                    w = 1
                    if rw_lb == REWEIGHT.NEIGHBOR_LABELS_HIST:
                        diff = np.linalg.norm(nb[i]-nb[j])
                        # tmp = deg[i] + deg[j] + deg[i]*deg[j]
                        tmp = 2
                        w = np.exp(-np.square(diff) / tmp)
                        #w = (deg[i] + deg[j]) / (1.0 + diff)
                    rs.append('{}\t{}\t{}\n'.format(i, j, w))
            ne = len(rs)
            if ne > 0 :
                # save to file
                N = len(gmap[gid])
                lb = glabel[gid]
                lb = lbs.index(lb)
                outfile = os.path.join(outfolder, 'nodes_{}_edges_{}_gid_{}_lb_{}_index_1_adj.nse'.format(N, ne, gid, lb))
                with open(outfile, 'w') as wf:
                    wf.writelines(rs)
            count += 1
            # if count > 2:
            #     break
        print('Saved {} files'.format(count))
    else:
        print('Unknow output format={} (should be .mat or .nse)'.format(ftype))