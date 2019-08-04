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
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.pipeline        import make_pipeline
from sklearn.preprocessing   import LabelEncoder
from sklearn.preprocessing   import StandardScaler
from sklearn.svm             import SVC

from collections import defaultdict
import operator

# for multiple kernel learning
import mklaren
import align
from mklaren.mkl.align import Align
from mklaren.mkl.alignf import Alignf

PARENT_PATH = r'F:\Research\ScaleVariant'
EXP_NAME    = 'exp_20190516'

MUTAG = 'MUTAG'
PFK_KERNEL = 4

"F:\Research\ScaleVariant\exp_20190516\MUTAG\gkernel\MUTAG_WL-VertexHist.txt"
"F:\Research\ScaleVariant\exp_20190516\MUTAG\gkernel\MUTAG_GraphletSampling.txt"
"F:\Research\ScaleVariant\exp_20190516\MUTAG\gkernel\MUTAG_ShortestPath.txt"
"F:\Research\ScaleVariant\exp_20190516\MUTAG\gkernel\MUTAG_VertexHist.txt"
"F:\Research\ScaleVariant\exp_20190516\MUTAG\gkernel\MUTAG_WL-ShortestPath.txt"
"F:\Research\ScaleVariant\exp_20190516\MUTAG\gkernel\MUTAG_WL-Subtree.txt"

FET_SCALE  = 'scale-variant'
FET_SUB_SCALE  = 'sub-scale-variant'
FET_GRAPHLET = 'GraphletSampling'
FET_SHORTEST = 'ShortestPath'
FET_WL_SH = 'WL-ShortestPath'
FET_WL_SB = 'WL-Subtree'
FET_WL_VH = 'WL-VertexHist'
FET_VH    = 'VertexHist'
FET_PROP_ATTR = 'PropagationAttr'
FET_GEOM  = 'GeometricRandomWalk'
FET_EXP   = 'ExponentialRandomWalk'
FET_THETA = 'SvmTheta'
FET_HB_LB = 'Hilbert_label_p1_2_p2_2'
FET_HB_ATT = 'Hilbert_attr_p1_2_p2_2'
FET_HB_LB_ATT = 'Hilbert_label_attr_p1_2_p2_2'


ker_methods = [FET_SCALE, FET_SUB_SCALE, FET_GRAPHLET, FET_SHORTEST, FET_WL_SB, FET_WL_SH, FET_WL_VH, FET_VH, FET_PROP_ATTR,\
    FET_GEOM, FET_EXP, FET_THETA, FET_HB_LB, FET_HB_ATT, FET_HB_LB_ATT]

def get_exp_path(parent_path, data_name, exp_name):
    exp_path = r'{}/{}/{}'.format(parent_path, exp_name, data_name)
    return exp_path

class KernelParams:
    def __init__(self, method, T1, T2, Tmax, thres, infval, t):
        self.method = method
        self.T1 = T1
        self.T2 = T2
        self.Tmax = Tmax
        self.thres = thres
        self.infval = infval
        self.t = t # for Fisher Persistence kernel

def normalize_kernel(kermat):
    #kermat[kermat <= 0] = 0
    # sm = np.sum(kermat.diagonal() <= 0)
    # if sm > 0:
    #     return None
    for i in range(kermat.shape[0]):
        if kermat[i, i] == 0:
            print('Kernel ill defined i =', i)
            kermat[i, i] = 1

    D = np.diag(1.0/np.sqrt(np.diag(kermat)))
    kermat = np.dot(np.dot(D, kermat), D)
    print('Normalized matrix')
    return kermat

def centered_kernel(kermat):
    '''
    Centered kernel aligmment
    '''
    if kermat is None:
        return None
    nz = kermat.shape[0]
    onev = np.ones(nz)
    onesub = (np.identity(nz) - np.dot(onev, np.transpose(onev)) / nz)
    kermat = np.dot(np.dot(onesub, kermat), onesub)
    return kermat


def get_kernel_path(folder, mt, exp_path, ph_name, dim, kprm):
    method, T1, T2, thres, infval, Tmax = kprm.method, kprm.T1, kprm.T2, kprm.thres, kprm.infval, kprm.Tmax
    if mt == FET_SCALE:
        kerpath = os.path.join(exp_path, '{}/{}_d_{}_method_{}_T1_{}_T2_{}_tmax_{}_thres_{}_inf_{}.txt'.format(folder, ph_name, \
            dim, method, T1, T2, Tmax, thres, infval))
    else:
        kerpath = 'DumpppppppNotFound'
    if not os.path.isfile(kerpath):
        print('Not found {}'.format(kerpath))
    return kerpath

def load_graph_kernel(data_name, exp_path, method, normalize):
    
    kerfile = os.path.join(exp_path, 'gkernel/{}_{}.txt'.format(data_name, method))
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
    print('Load kernel', kermat.shape)
    return kermat

def load_data(filename):
    """
    Read from setting file
    """
    lbs = []
    datls = []
    
    with open(filename, 'r', encoding='utf-8') as rf:
        lines = rf.readlines()
        print(len(lines))
        for line in lines:
            lb = int(re.search('lb_([0-9]+)_', line).groups()[0])
            lbs.append(lb)
            datls.append(len(lbs)-1)
        rf.close()

    X_index = np.array(datls).astype(np.int32)
    y = np.array(lbs).astype(np.int32)
    print('Loaded setting file, shape ', X_index.shape, y.shape)
    return X_index, y

def svc_classify(X, y, train_index, test_index, mt, opt_C=None):
    if mt == 'common':
        X_train, X_test = X[np.ix_(train_index)], X[np.ix_(test_index)]
    else:
        X_train, X_test = X[np.ix_(train_index, train_index)], X[np.ix_(test_index, train_index)]
    y_train, y_test = y[np.ix_(train_index)], y[np.ix_(test_index)]

    #C_grid = [1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0, 2e0, 5e0, 1e1, 2e1, 5e1, 1e2]
    C_grid = [1e-2, 1e-1, 1e0, 1e1, 1e2]
    #C_grid = (10. ** np.arange(1,10,1) / len(y_train)).tolist()
    if mt == 'common':
        best_clf = GridSearchCV(SVC(), cv=5, param_grid={'kernel':('linear', 'rbf'), 'C': C_grid})
    else:
        if opt_C is None:
            best_clf = GridSearchCV(SVC(kernel='precomputed'), cv=5, param_grid={"C": C_grid})
        else:
            best_clf = SVC(kernel='precomputed', C=opt_C)
    
    best_clf.fit(X_train, y_train)

    train_sc = best_clf.score(X_train, y_train)
    test_sc  = best_clf.score(X_test, y_test)
    
    # Compute confusion matrix
    # y_pred = best_clf.predict(X_test)
    # cnf_matrix = confusion_matrix(y_test, y_pred)
    return train_sc, test_sc

def add_kernel(K1, K2):
    if K1 is None:
        return K2
    if K2 is None:
        return K1
    return (K1 + K2)

def multiply_kernel(K1, K2):
    if K1 is None:
        return K2
    if K2 is None:
        return K1
    return np.multiply(K1, K2)

def combine_kernel(K1, K2, alpha):
    if K1 is None:
        return K2
    if K2 is None:
        return K1
    kX = alpha*K1 + (1.0-alpha)*K2
    return kX

def select_kernel(y, train_index, K1, K2):
    skf = StratifiedKFold(n_splits=5, shuffle=False, random_state=0)
    alpha_grid = [0.0, 1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 8e-1, 9e-1, 95e-2, 98e-2, 99e-2, 999e-3, 1.0]
    #alpha_grid = [1.0]
    C_grid = [1e-2, 1e-1, 1e0, 1e1, 1e2]
    #C_grid = [1.0]

    scores = defaultdict(np.float32)
    for alpha in alpha_grid:
        kX = combine_kernel(K1, K2, alpha)
        if kX is None:
            return None, 0.0
        X_train = kX[np.ix_(train_index, train_index)]
        y_train = y[np.ix_(train_index)]

        for C in C_grid:
            local_test_sc = []
            for local_train_idx, local_test_idx in skf.split(X_train, y_train):
                #print(local_train_idx, local_test_idx)
                X_local_train = X_train[np.ix_(local_train_idx, local_train_idx)]
                X_local_test  = X_train[np.ix_(local_test_idx, local_train_idx)]
                y_local_train, y_local_test = y[np.ix_(local_train_idx)], y[np.ix_(local_test_idx)]
                local_clf = SVC(kernel='precomputed', C=C)
                local_clf.fit(X_local_train, y_local_train)

                local_test_sc.append(local_clf.score(X_local_test, y_local_test))

            scores[(alpha, C)] = np.mean(np.array(local_test_sc))
            #print('alpha={}, C={}, scores={}'.format(alpha, C, scores[(alpha, C)] ))
        if (K1 is None) or (K2 is None):
            print('Without combination')
            break
    #mkey = max(scores.items(), key=operator.itemgetter(1))[0]
    mkey = max(scores, key=scores.get)
    opt_alpha, opt_C = mkey[0], mkey[1]
    print('Opt (alpha, C)', mkey, 'with score={}'.format(scores[mkey]))
    kX = combine_kernel(K1, K2, opt_alpha)
    return kX, opt_C

def multiple_kernel_learning(y, train_index, K1, K2, label="normal"):
    if K1 is None:
        return K2
    if K2 is None:
        return K1
    K1_train = K1[np.ix_(train_index, train_index)]
    K2_train = K2[np.ix_(train_index, train_index)]
    
    y_train = y[np.ix_(train_index)]

    model = Alignf(typ="convex")
    #model = Align()
    model.fit([K1_train, K2_train], y_train)
    mu = model.mu
    #combined_kernel = lambda x, y: \
    #    mu[0] * K1(x, y) + mu[1] * K2(x, y)
    print(label, mu)
    #combine_kernel = mu[0] * centered_kernel(K1) + mu[1] * centered_kernel(K2)
    combine_kernel = mu[0] * K1 + mu[1] * K2
    return combine_kernel

np.set_printoptions(precision=5)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataname', '-d', type=str, default=MUTAG)
    parser.add_argument('--log', type=str, default='log')
    parser.add_argument('--parentpath', type=str, default=PARENT_PATH)
    parser.add_argument('--expname', type=str, default=EXP_NAME)
    parser.add_argument('--phname', type=str, default='ph_20180628_norm_0')
    parser.add_argument('--sphname', type=str, default='ph_20180628_norm_0')
    parser.add_argument('--method', '-me', type=int, default=0)
    parser.add_argument('--T1', type=float, default=0.0)
    parser.add_argument('--T2', type=float, default=1.0)
    parser.add_argument('--Tmax', type=float, default=0.0)
    parser.add_argument('--thres0', type=float, default=0.0)
    parser.add_argument('--thres1', type=float, default=0.0)
    parser.add_argument('--sTmax', type=float, default=0.0)
    parser.add_argument('--sthres0', type=float, default=0.0)
    parser.add_argument('--sthres1', type=float, default=0.0)
    parser.add_argument('--weight', '-w', type=float, default=0.5)
    parser.add_argument('--time', type=float, default=1.0)
    parser.add_argument('--infval', type=float, default=0.0)
    parser.add_argument('--norm', type=int, default=1)
    parser.add_argument('--nums', '-nu', type=int, default=1)
    parser.add_argument('--combine', '-cb', type=int, default=-1) # 0: add, 1: multiply, other: multiple learning
    parser.add_argument('--sp', type=int, default=0) # specific index of kernel method, -1 for test all methods
    parser.add_argument('--sb', type=int, default=-1) # sepecific index of subkernel
    parser.add_argument('--dp', type=int, default=-1) # specific dim for calculating kernel, -1 for use all dimensions
    parser.add_argument('--sdp', type=int, default=-1) # specific dim for calculating sub kernel, -1 for use all dimensions
    
    args = parser.parse_args()

    parent_path, data_name, exp_name, ph_name = args.parentpath, args.dataname, args.expname, args.phname
    sph_name, sdp = args.sphname, args.sdp
    dp, sp, sb, norm, infval, weight = args.dp, args.sp, args.sb, args.norm, args.infval, args.weight
    combine = args.combine

    dir_path = os.path.dirname(os.path.realpath(__file__))
    log_path = os.path.join(dir_path, args.log)
    exp_path = get_exp_path(parent_path, data_name, exp_name)

    if os.path.exists(log_path) == False:
        os.makedirs(log_path)

    log_filename = '{}_{}_{}_T1_{}_T2_{}_thres0_{}_thres1_{}_nums_{}_method_{}_norm_{}_infval_{}_t_{}_sp_{}_sub_{}_comb_{}.log'.format(
        data_name, exp_name, ph_name,
        args.T1, args.T2, args.thres0, args.thres1,
        args.nums, args.method, args.norm, args.infval, args.time, sp, sb, combine
    )
    log_filename = os.path.join(log_path, log_filename)
    logger = get_module_logger(__name__, log_filename)
    logger.info(log_filename)

    setting_file = os.path.join(exp_path, 'matlist_grakel.txt')
    X_index, y = load_data(setting_file)

    # Load kernels
    sp_mth, sb_mth = '', ''
    if sp >= 0 and sp < len(ker_methods):
        sp_mth = ker_methods[sp]
    if sb > 0 and sb < len(ker_methods):
        sb_mth = ker_methods[sb]
    kerX = defaultdict()

    # scale-variant type
    kX0, kX1 = None, None
    skX0, skX1 = None, None
    kX, skX = None, None

    for mt in ker_methods:
        if len(sp_mth) > 0 and mt != sp_mth:
            continue
        elif mt == FET_SCALE:
            kprm0 = KernelParams(args.method, args.T1, args.T2, args.Tmax, args.thres0, args.infval, args.time)
            kprm1 = KernelParams(args.method, args.T1, args.T2, args.Tmax, args.thres1, args.infval, args.time)

            kX0 = load_kernel('kernel', mt, exp_path, ph_name, dp, 0, kprm0, norm)
            kX1 = load_kernel('kernel', mt, exp_path, ph_name, dp, 1, kprm1, norm)
            
            kX = add_kernel(kX0, kX1)
            if kX is not None:
                print('Scale kernel loaded mt={} with shape'.format(mt), kX.shape)
            # load subkernel
            skX = None
            if sb_mth == FET_SUB_SCALE:
                kprm0 = KernelParams(args.method, args.T1, args.T2, args.sTmax, args.sthres0, args.infval, args.time)
                kprm1 = KernelParams(args.method, args.T1, args.T2, args.sTmax, args.sthres1, args.infval, args.time)
                skX0 = load_kernel('kernel', mt, exp_path, sph_name, sdp, 0, kprm0, norm)
                skX1 = load_kernel('kernel', mt, exp_path, sph_name, sdp, 1, kprm1, norm)
                skX = add_kernel(skX0, skX1)
            else:
                skX = load_graph_kernel(data_name, exp_path, sb_mth, norm)
            if skX is not None:
                print('Sub kernel loaded mt={} with shape'.format(sb_mth), skX.shape)
                kX  = add_kernel(kX, skX)
                #kX = np.multiply(kX, kX) + np.multiply(skX, skX)
            if kX is not None:
                kerX[mt] = kX
        elif mt in ker_methods:
            # load graph kernel
            gX = load_graph_kernel(data_name, exp_path, mt, norm)
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
                if mt == FET_SCALE:
                    kX = multiple_kernel_learning(y, train_index, kX0, kX1, label="learning main kernels")
                    if sb_mth == FET_SUB_SCALE:
                        skX = multiple_kernel_learning(y, train_index, skX0, skX1, label="learning sub kernels")
                    if combine == 0:
                        #print('Add kernel')
                        kerX[mt] = add_kernel(kX, skX)
                    elif combine == 1:
                        #print('Multiply kernel')
                        kerX[mt] = multiply_kernel(kX, skX)
                    else:
                        kerX[mt] = multiple_kernel_learning(y, train_index, kX, skX, label="learning combine kernels")
                    #kerX[mt] = normalize_kernel(kerX[mt])
                    #model = multiple_kernel_learning(y, train_index, kX0, kX1, skX0, skX1)
                
                train_sc, test_sc = svc_classify(kerX[mt], y, train_index, test_index, mt)
                local_train.append(train_sc)
                local_test.append(test_sc)
                #logger.debug('n={}, mt={}, score: train={}, test={}'.format(n, mt, train_sc, test_sc))

            avg_test_local  = np.mean(local_test)
            avg_train_local = np.mean(local_train)

            global_test[mt].append(avg_test_local)
            global_train[mt].append(avg_train_local)

            #logger.debug('n={}, mean local score: train={}, test={}'.format(n, avg_train_local, avg_test_local))
        
            avg_test_global,  std_test_global  = 100*np.mean(global_test[mt]),  100*np.std(global_test[mt])
            avg_train_global, std_train_global = 100*np.mean(global_train[mt]), 100*np.std(global_train[mt])

            logger.debug('n={}, mt={}, glocal score (mean, std): train=( {}, {} ), test=( {}, {} )'.format(n, mt,\
                avg_train_global, std_train_global,
                avg_test_global, std_test_global))
        logger.debug('')
            