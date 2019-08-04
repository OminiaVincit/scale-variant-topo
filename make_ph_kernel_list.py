import numpy as np
import glob
import os
import argparse
import re
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exppath', '-p', type=str, required=True)
    parser.add_argument('--infolder', '-i', type=str, required=True)
    parser.add_argument('--dim', '-d', type=int, default=1)
    args = parser.parse_args()
    print(args)

    parent = args.exppath
    infolder, dim = args.infolder, args.dim

    outfile = r'{}\barlist_{}_d_{}.txt'.format(parent, infolder, dim)
    lbfile   = r'{}\labels_{}_d_{}.txt'.format(parent, infolder, dim)

    filedict = defaultdict()

    for phfile in glob.glob(r"{}\{}\diffusion_barcode_nodes_*_dim_{}.txt".format(parent, infolder, dim)):
        lb = int(re.search('lb_([0-9]+)_', phfile).groups()[0])
        gid = int(re.search('gid_([0-9]+)_', phfile).groups()[0])
        filedict[gid] = phfile

    filels = []
    for k, v in sorted(filedict.items()):
        filels.append(v + '\n')
    
    print('Num files={}'.format(len(filels)))
    with open(outfile, 'w') as wf:
        wf.writelines(filels)
        wf.close()
    
    # with open(lbfile, 'w') as wf:
    #     wf.writelines(lbs)
    #     wf.close()