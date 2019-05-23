import numpy as np
import glob
import os
import argparse
import re
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exppath', '-e', type=str, required=True)
    parser.add_argument('--infolder', '-i', type=str, default='grakel')
    args = parser.parse_args()
    print(args)

    parent, infolder = args.exppath, args.infolder

    pathlist = [os.path.join(parent, o) for o in os.listdir(parent) if os.path.isdir(os.path.join(parent, o))]

    for p in pathlist:
        outfile = r'{}\matlist_{}.txt'.format(p, infolder)
        filedict = defaultdict()
        matpath = os.path.join(p, infolder)
        if os.path.isdir(matpath) == False:
            continue

        for matfile in glob.glob(r'{}\nodes_*.mat'.format(matpath)):
            lb = int(re.search('lb_([0-9]+)_', matfile).groups()[0])
            gid = int(re.search('gid_([0-9]+)_', matfile).groups()[0])
            filedict[gid] = matfile

        filels = []
        for k, v in sorted(filedict.items()):
            filels.append(v + '\n')
        
        print('{} Num files={}'.format(p, len(filels)))
        with open(outfile, 'w') as wf:
            wf.writelines(filels)
            wf.close()