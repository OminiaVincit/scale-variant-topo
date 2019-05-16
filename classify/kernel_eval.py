import numpy as np
from loginit import get_module_logger

from sklearn.decomposition import PCA, KernelPCA
from sklearn.model_selection import GridSearchCV, StratifiedShuffleSplit, StratifiedKFold
from sklearn.model_selection import cross_validate
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import confusion_matrix, roc_auc_score
import re
import sys
import os
import time
import argparse
import random
import itertools
import cmocean

import warnings
import sklearn.exceptions
#warnings.filterwarnings("ignore", category=sklearn.exceptions.UndefinedMetricWarning)

import matplotlib as mpl
import matplotlib.pyplot as plt

from sklearn.base            import clone
from sklearn.metrics         import accuracy_score
from sklearn.metrics         import average_precision_score
from sklearn.metrics         import classification_report
from sklearn.metrics         import f1_score
from sklearn.metrics         import precision_score
from sklearn.metrics         import recall_score
from sklearn.metrics         import roc_auc_score
from sklearn.model_selection import ParameterGrid
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.pipeline        import make_pipeline
from sklearn.preprocessing   import LabelEncoder
from sklearn.preprocessing   import StandardScaler
from sklearn.svm             import SVC

from collections import defaultdict

PARENT_PATH = r'F:\Research\ScaleVariant'
EXP_NAME    = 'exp_20190516'

MUTAG = 'MUTAG'
PFK_KERNEL = 4

FET_SCALE  = 'scale-variant'
FET_COMMON = 'common'
FET_GRAPHLET = 'Graphlet'
FET_KSTEP = 'KStepRandomWalk'
FET_GEOM  = 'GeometricRandomWalk'
FET_EXP   = 'ExponentialRandomWalk'
FET_SHORTEST = 'ShortestPath'
FET_WL = 'WL'
FET_AVG_SCALE = 'scale_avg'
FET_AVG_SCALE_NORM = 'scale_norm_avg'

ker_methods = [FET_SCALE,FET_GRAPHLET, FET_KSTEP, FET_GEOM, FET_EXP, FET_SHORTEST, FET_WL, FET_COMMON, FET_AVG_SCALE, FET_AVG_SCALE_NORM]
ker_pal = {
    FET_GRAPHLET: '4', 
    FET_KSTEP : '2', 
    FET_GEOM : '0.05', 
    FET_EXP : '0.1', 
    FET_SHORTEST : 'NA',
    FET_WL : '5'
}

def get_exp_path(parent_path, data_name, exp_name):
    exp_path = r'{}/{}/{}'.format(parent_path, exp_name, data_name)
    return exp_path

class KernelParams:
    def __init__(self, method, T1, T2, thres, infval, t):
        self.method = method
        self.T1 = T1
        self.T2 = T2
        self.thres = thres
        self.infval = infval
        self.t = t # for Fisher Persistence kernel

def normalize_kernel(kermat):
    #kermat[kermat <= 0] = 0
    sm = np.sum(kermat.diagonal() <= 0)
    if sm > 0:
        return None

    D = np.diag(1.0/np.sqrt(np.diag(kermat)))
    kermat = np.dot(np.dot(D, kermat), D)
    print('Normalized matrix')
    return kermat

def get_kernel_path(folder, mt, exp_path, ph_name, dim, kprm):
    method, T1, T2, thres, infval = kprm.method, kprm.T1, kprm.T2, kprm.thres, kprm.infval
    if mt == FET_SCALE:
        kerpath = os.path.join(exp_path, '{}/{}_d_{}_method_{}_T1_{}_T2_{}_tmax_{}_thres_{}_inf_{}.txt'.format(folder, ph_name, \
            dim, method, T1, T2, 0.0, thres, infval))
    elif mt == FET_AVG_SCALE or mt == FET_AVG_SCALE_NORM:
        kerpath = os.path.join(exp_path, '{}/{}_method_{}_d_{}.txt'.format(folder, mt, method, dim))
    else:
        kerpath = 'DumpppppppNotFound'
    if not os.path.isfile(kerpath):
        print('Not found {}'.format(kerpath))
    return kerpath

def load_graph_kernel(folder, exp_path, method, par, normalize):
    kerfile = os.path.join(exp_path, '{}/graph_kernel_{}_par_{}.txt'.format(folder, method, par))
    if method == FET_GRAPHLET and os.path.isfile(kerfile) == False:
        kerfile = os.path.join(exp_path, '{}/graph_kernel_{}_par_3.txt'.format(folder, method))
    if os.path.isfile(kerfile) == False:
        print('Not found graph kernel', method)
        return None
    kermat = np.loadtxt(kerfile)
    if kermat is None:
        return None
    if normalize > 0:
        kermat = normalize_kernel(kermat)
    return kermat

def load_kernel(folder, mt, exp_path, ph_name, dp, dim, kprm, normalize):
    if dp >= 0 and dp != dim:
        print ('Specified dim not match dp={}, dim={}'.format(dp, dim))
        return
    kerfile = get_kernel_path(folder, mt, exp_path, ph_name, dim, kprm)
    if os.path.isfile(kerfile) == False:
        return None
    kermat = np.loadtxt(kerfile, dtype=np.float32)
    if kermat is None:
        return None
    print('Dim=', dim, 'Kermat shape', kermat.shape)
    if kprm.method == PFK_KERNEL:
        kermat = kermat / float(-kprm.t)
        kermat = np.exp(kermat)
    if normalize > 0:
        kermat = normalize_kernel(kermat)
    return kermat

def load_data(filename):
    """
    Read from setting file
    """
    lbs = []
    datls = []
    
    with open(filename, 'r') as rf:
        lines = rf.readlines()
        for line in lines:
            lb = int(re.search('lb_([0-9]+)_', line).groups()[0])
            lbs.append(lb)
            datls.append(len(lbs)-1)
        rf.close()

    X_index = np.array(datls).astype(np.int32)
    y = np.array(lbs).astype(np.int32)

    return X_index, y

def svc_classify(X, y, train_index, test_index, mt):
    if mt == 'common':
        X_train, X_test = X[np.ix_(train_index)], X[np.ix_(test_index)]
    else:
        X_train, X_test = X[np.ix_(train_index, train_index)], X[np.ix_(test_index, train_index)]
    y_train, y_test = y[np.ix_(train_index)], y[np.ix_(test_index)]

    if mt == 'common':
        best_clf = GridSearchCV(SVC(), cv=5, param_grid={'kernel':('linear', 'rbf'), 'C':[1e-2, 1e-1, 1e0, 1e1, 1e2]})
    else:
        best_clf = GridSearchCV(SVC(kernel='precomputed'), cv=5, param_grid={"C": [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]})
    best_clf.fit(X_train, y_train)

    train_sc = best_clf.score(X_train, y_train)
    test_sc  = best_clf.score(X_test, y_test)
    
    # Compute confusion matrix
    # y_pred = best_clf.predict(X_test)
    # cnf_matrix = confusion_matrix(y_test, y_pred)
    return train_sc, test_sc

np.set_printoptions(precision=5)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname', '-d', type=str, default=MUTAG)
    parser.add_argument('--log', type=str, default='log')
    parser.add_argument('--parentpath', type=str, default=PARENT_PATH)
    parser.add_argument('--expname', type=str, default=EXP_NAME)
    parser.add_argument('--phname', type=str, default='ph_20180628_norm_0')
    parser.add_argument('--method', '-me', type=int, default=0)
    parser.add_argument('--T1', type=float, default=0.0)
    parser.add_argument('--T2', type=float, default=1.0)
    parser.add_argument('--thres0', type=float, default=0.0)
    parser.add_argument('--thres1', type=float, default=0.0)
    parser.add_argument('--time', type=float, default=1.0)
    parser.add_argument('--infval', type=float, default=0.0)
    parser.add_argument('--norm', type=int, default=1)
    parser.add_argument('--nums', '-nu', type=int, default=1)
    parser.add_argument('--sp', type=int, default=0) # specific index of kernel method, -1 for test all methods
    parser.add_argument('--dp', type=int, default=-1) # specific dim for calculating kernel, -1 for use all dimensions
    args = parser.parse_args()

    parent_path, data_name, exp_name, ph_name = args.parentpath, args.dataname, args.expname, args.phname
    dp, sp, norm, infval = args.dp, args.sp, args.norm, args.infval

    dir_path = os.path.dirname(os.path.realpath(__file__))
    log_path = os.path.join(dir_path, args.log)
    exp_path = get_exp_path(parent_path, data_name, exp_name)

    if os.path.exists(log_path) == False:
        os.makedirs(log_path)

    log_filename = '{}_{}_{}_T1_{}_T2_{}_thres0_{}_thres1_{}_nums_{}_method_{}_norm_{}_infval_{}_t_{}.log'.format(
        data_name, exp_name, ph_name,
        args.T1, args.T2, args.thres0, args.thres1,
        args.nums, args.method, args.norm, args.infval, args.time
    )
    log_filename = os.path.join(log_path, log_filename)
    logger = get_module_logger(__name__, log_filename)
    logger.info(log_filename)

    # Load kernels
    sp_mth = ''
    if sp >= 0 and sp < len(ker_methods):
        sp_mth = ker_methods[sp]

    kerX = defaultdict()

    for mt in ker_methods:
        if len(sp_mth) > 0 and mt != sp_mth:
            continue
        if mt == FET_COMMON:
            # features type
            feature_file = os.path.join(exp_path, 'kernel/common_features.npy')
            if os.path.isfile(feature_file):
                fetX = np.load(feature_file)
                fetX[np.isnan(fetX)] = 0.0
                if norm > 0:
                    for i in range(fetX.shape[1]):
                        v = fetX[:, i]
                        vmin, vmax = v.min(), v.max()
                        if vmin < vmax:
                            fetX[:, i] = (v - vmin) / (vmax - vmin)
                kerX[mt] = fetX
                print('Features loaded')
        elif mt == FET_SCALE or mt == FET_AVG_SCALE or mt == FET_AVG_SCALE_NORM:
            kprm0 = KernelParams(args.method, args.T1, args.T2, args.thres0, args.infval, args.time)
            kprm1 = KernelParams(args.method, args.T1, args.T2, args.thres1, args.infval, args.time)

            # scale-variant type
            kX0, kX1 = None, None

            kX0 = load_kernel('kernel', mt, exp_path, ph_name, dp, 0, kprm0, norm)
            kX1 = load_kernel('kernel', mt, exp_path, ph_name, dp, 1, kprm1, norm)

            setting_file = os.path.join(exp_path, 'barlist_{}_d_0.txt'.format(ph_name))

            X_index, y = load_data(setting_file)
            kX = kX0
            if kX1 is not None:
                if kX is not None:
                    kX = kX + kX1
                    print('Combined kernel')
                else:
                    kX = kX1
            if kX is not None:
                kerX[mt] = kX
                print('Scale kernel loaded mt={}'.format(mt))
                print('Shape = ', kX.shape)
        elif mt in ker_pal:
            # load graph kernel
            gX = load_graph_kernel('kernel', exp_path, mt, ker_pal[mt], norm)
            if gX is not None:
                kerX[mt] = gX
                print('{} kernel loaded'.format(mt))

    np.set_printoptions(precision=5)

    
    global_test, global_train = defaultdict(), defaultdict()
    for mt in kerX.keys():
        global_test[mt] = []
        global_train[mt]  = []
    
    # repeat for nums
    for n in range(args.nums):
        # ten-fold
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=n)
        for mt in kerX.keys():
            local_test, local_train = [], []
            for train_index, test_index in skf.split(X_index, y):
                train_sc, test_sc = svc_classify(kerX[mt], y, train_index, test_index, mt)
                local_train.append(train_sc)
                local_test.append(test_sc)
                #logger.debug('n={}, mt={}, score: train={}, test={}'.format(n, mt, train_sc, test_sc))

            avg_test_local  = np.mean(local_test)
            avg_train_local = np.mean(local_train)

            global_test[mt].append(avg_test_local)
            global_train[mt].append(avg_train_local)

            #logger.debug('n={}, mean local score: train={}, test={}'.format(n, avg_train_local, avg_test_local))
        
            avg_test_global,  std_test_global  = np.mean(global_test[mt]),  np.std(global_test[mt])
            avg_train_global, std_train_global = np.mean(global_train[mt]), np.std(global_train[mt])

            logger.debug('n={}, mt={}, glocal score (mean, std): train=({},{}), test=({},{})'.format(n, mt,\
                avg_train_global, std_train_global,
                avg_test_global, std_test_global))
        logger.debug('')
            