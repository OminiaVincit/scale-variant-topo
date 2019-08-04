import numpy as np
import glob
import os
import argparse
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--parent', '-p', type=str, required=True)
    parser.add_argument('--data', '-d', type=str, default='BA')
    parser.add_argument('--midx', '-m', type=int, default=10)
    args = parser.parse_args()
    print(args)
    parent, m = args.parent, args.midx

    dks = ['original', 'dk2.5', 'dk2.1', 'dk2.0', 'dk1.0']

    for dim in [0, 1]:
        outfile = r'{}\barlist_join_{}_d_{}.txt'.format(parent, m, dim)
        filels = []
        for dk in dks:
            infolder = os.path.join(parent, 'ph_{}_{}'.format(args.data, dk))
            print('Get files from ', infolder)
            for phfile in glob.glob(r"{}\diffusion_barcode_*_dim_{}.txt".format(infolder, dim)):
                index = int(re.search('index_([0-9]+)_', phfile).groups()[0])
                if index <= m:
                    filels.append(phfile + '\n')
                    
        print('dim={}, Num files={}'.format(dim, len(filels)))
        with open(outfile, 'w') as wf:
            wf.writelines(filels)
            wf.close()
    