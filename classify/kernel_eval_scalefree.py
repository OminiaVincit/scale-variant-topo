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


EXPPATH = r'F:\Research\ScaleVariant'

class KernelParams:
    def __init__(self, method, T1, T2, thres, infval, interval, Tmax):
        self.method = method
        self.T1 = T1
        self.T2 = T2
        self.thres = thres
        self.infval = infval
        self.interval = interval
        self.Tmax = Tmax

def normalize_kernel(kermat):
    #kermat[kermat <= 0] = 0
    sm = np.sum(kermat.diagonal() <= 0)
    if sm > 0:
        return None

    D = np.diag(1.0/np.sqrt(np.diag(kermat)))
    kermat = np.dot(np.dot(D, kermat), D)
    print('Normalized matrix')
    return kermat

def get_kpca_gen_file(folder, exp_path, setting_file, kprm):
    T1, T2, tmax, interval, thres, infval = kprm.T1, kprm.T2, kprm.Tmax, kprm.interval, kprm.thres, kprm.infval
    kername = setting_file.replace('barlist_', '')
    kername = kername.replace('.txt', '_method_0_T1_{}_T2_{}_tmax_{}_interval_{}_thres_{}_inf_{}.txt'.format(T1, T2, tmax, interval, thres, infval))
    kerpath = os.path.join(exp_path, folder)
    kerpath = os.path.join(kerpath, kername)
    return kerpath


def load_kernel(folder, exp_path, setting_file, kprm, normalize):
    kerfile = get_kpca_gen_file(folder, exp_path, setting_file, kprm)
    if os.path.isfile(kerfile) == False:
        print('Not found', kerfile)
        return None
    kermat = np.loadtxt(kerfile, dtype=np.float32)
    if kermat is None:
        return None
    for i in range(kermat.shape[0]):
        if kermat[i, i] == 0:
            print('Kernel ill defined i =', i)
            kermat[i, i] = 1
    if normalize > 0:
        kermat = normalize_kernel(kermat)
    return kermat

def load_data(filepath):
    """
    Reading setting file
    """
    lbs = []
    datls = []
    with open(filepath, 'r') as rf:
        lines = rf.readlines()
        for line in lines:
            lb = -1
            if 'BA' in line:
                lb = 1
            elif 'conf' in line:
                lb = 2
            if lb > 0:
                lbs.append(lb)
                datls.append(len(lbs)-1)
        rf.close()
    print('Number of files = {}'.format(len(lbs)))
    X_index = np.array(datls).astype(np.int32)
    y = np.array(lbs).astype(np.int32)

    return X_index, y

def svc_classify(X, y, train_index, test_index):
    X_train, X_test = X[np.ix_(train_index, train_index)], X[np.ix_(test_index, train_index)]
    y_train, y_test = y[np.ix_(train_index)], y[np.ix_(test_index)]

    #C_grid = [1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0, 2e0, 5e0, 1e1, 2e1, 5e1, 1e2]
    C_grid = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
    #C_grid = (10. ** np.arange(1,10,1) / len(y_train)).tolist()
    best_clf = GridSearchCV(SVC(kernel='precomputed'), cv=5, param_grid={"C": C_grid})
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
    parser.add_argument('--log', type=str, default='logscale')
    parser.add_argument('--exppath', type=str, default=EXPPATH)
    parser.add_argument('--setting', type=str, default='barlist_scalefree_d_1.txt')
    parser.add_argument('--method', '-me', type=int, default=0)
    parser.add_argument('--T1', type=float, default=0.0)
    parser.add_argument('--T2', type=float, default=1.0)
    parser.add_argument('--Tmax', type=float, default=100.0)
    parser.add_argument('--thres', type=float, default=0.0)
    parser.add_argument('--interval', type=int, default=1)
    parser.add_argument('--infval', type=float, default=0.0)
    parser.add_argument('--norm', type=int, default=1)
    parser.add_argument('--nums', '-nu', type=int, default=1)
    args = parser.parse_args()

    exppath, settingname = args.exppath, args.setting
    T1, T2, Tmax, thres, interval, infval, norm, nums = args.T1, args.T2, args.Tmax, args.thres, args.interval, args.infval, args.norm, args.nums
    dir_path = os.path.dirname(os.path.realpath(__file__))
    log_path = os.path.join(dir_path, args.log)

    if os.path.exists(log_path) == False:
        os.makedirs(log_path)

    log_filename = '{}_T1_{}_T2_{}_Tmax_{}_thres_{}_interval_{}_nums_{}_norm_{}_infval_{}.log'.format(
        settingname.replace('.txt', ''), T1, T2, Tmax, thres, interval, nums, norm, infval
    )
    log_filename = os.path.join(log_path, log_filename)
    logger = get_module_logger(__name__, log_filename)
    logger.info(log_filename)

    setting_file = os.path.join(exppath, settingname)
    selected_lbs = None

    # Load kernels
    kprm = KernelParams(0, T1, T2, thres, infval, interval, Tmax)
    kerX = load_kernel('kernel', exppath, settingname, kprm, norm)
    print('kerX shape', kerX.shape)
    global_test, global_train = [], []

    X_index, y = load_data(setting_file)

    # repeat for nums
    for n in range(args.nums):
        local_test, local_train = [], []
        # 2-fold
        skf = StratifiedKFold(n_splits=2, shuffle=True, random_state=n)
        for train_index, test_index in skf.split(X_index, y):
            train_sc, test_sc = svc_classify(kerX, y, train_index, test_index)  
            local_train.append(train_sc)
            local_test.append(test_sc)
        
        avg_test_local  = np.mean(local_test)
        avg_train_local = np.mean(local_train)

        global_train.append(avg_train_local)
        global_test.append(avg_test_local)

        logger.debug('n={}, local score: train={}, test={}'.format(n, avg_train_local, avg_test_local))
    
        avg_test_global,  std_test_global  = 100*np.mean(global_test),  100*np.std(global_test, dtype=np.float32)
        avg_train_global, std_train_global = 100*np.mean(global_train), 100*np.std(global_train, dtype=np.float32)

        logger.debug('n={}, glocal score (mean, std): train=( {}, {} ), test=( {}, {} )'.format(n, \
            avg_train_global, std_train_global,
            avg_test_global, std_test_global))
        logger.debug('')
    logger.debug('Final,{},{}'.format(avg_test_global, std_test_global)) 
