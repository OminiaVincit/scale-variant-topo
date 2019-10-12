import re
import numpy as np
import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--logfolder', type=str, default='logdkconf')
    parser.add_argument('--std', type=int, default=0)
    args = parser.parse_args()

    # tlist = [1.0]
    # for x in range(20):
    #     tlist.append(5*x+5)

    tlist = list(range(20, 51))
    for Tmax in tlist:
        rs = []
        for interval in [1]:
            filename = '{}/barlist_config_join_10_d_1_T1_0.0_T2_1.0_Tmax_{}_thres_0.0_interval_{}_nums_100_norm_1_infval_0.0.log'.format(
                args.logfolder, format(Tmax, '.1f'), interval)
            tmp = 'a'
            if os.path.isfile(filename) == True:
                with open(filename, 'r') as rf:
                    lines = rf.readlines()
                    for line in lines:
                        line = line.strip()
                        if 'Final' in line:
                            ls = line.split(',')
                            mstd = ls[-1]
                            acc = ls[-2]
                            if args.std > 0:
                                tmp = mstd
                            else:
                                tmp = acc
            rs.append(tmp)
        print(Tmax, ','.join(rs))