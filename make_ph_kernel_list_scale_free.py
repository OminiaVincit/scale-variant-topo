import numpy as np
import glob
import os
import argparse
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--parent', '-p', type=str, required=True)
    parser.add_argument('--data', '-d', type=str, default='scalefree')
    args = parser.parse_args()
    print(args)
    parent= args.parent

    dks = ['BA', 'conf']

    for dim in [0, 1]:
        outfile = r'{}\barlist_scalefree_d_{}.txt'.format(parent, dim)
        filels = []
        for dk in dks:
            infolder = os.path.join(parent, 'ph_{}_{}'.format(args.data, dk))
            print('Get files from ', infolder)
            for phfile in glob.glob(r"{}\diffusion_barcode_*_dim_{}.txt".format(infolder, dim)):
                filels.append(phfile + '\n')
                    
        print('dim={}, Num files={}'.format(dim, len(filels)))
        with open(outfile, 'w') as wf:
            wf.writelines(filels)
            wf.close()
    