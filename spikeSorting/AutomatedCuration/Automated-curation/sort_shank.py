import numpy as np
import joblib
from fastai.text.all import *
from prepare_features_MC import prepare_features_MC
from uptade_all import update_all2




def predict2(feat_mat):
    #clf = load_learner('incep_25_10.pkl')
    import pathlib
    posix_backup = pathlib.PosixPath
    try:
        pathlib.PosixPath = pathlib.WindowsPath
        clf = load_learner('C:\\Users\\Cornell\\Automated-curation\\MC_v1.pkl')
    finally:
        pathlib.PosixPath = posix_backup   
    
    feat_mat = tensor(feat_mat)
    feat_mat[:, 4, 0:83] += -0.5
    dl = clf.dls.test_dl(feat_mat)
    preds = clf.get_preds(dl=dl)
    a = preds[0].numpy()
    return a

def load_data(filebase, shank_num):
    f = open(filebase+'/info')
    info = f.read()
    lst_info = info.split('\n')
    session_name = lst_info[0]

    clu = np.load(filebase + '/' + session_name + '.clu.' + shank_num + '.npy')
    res = np.load(filebase + '/' + session_name + '.res.' + shank_num + '.npy')
    cc = np.load(filebase + '/' + session_name + '.cc.' + shank_num + '.npy')
    mspk = np.load(filebase + '/' + session_name + '.mspk.' + shank_num + '.npy')
    sspk = np.load(filebase + '/' + session_name + '.sspk.' + shank_num + '.npy')
    nspk_vec = np.load(filebase + '/' + session_name + '.nspk_vec.' + shank_num + '.npy')
    time_mat = np.load(filebase + '/' + session_name + '.timeMat.' + shank_num + '.npy')
    nspk_vec = np.squeeze(nspk_vec)
    return clu, res, cc, mspk, sspk, nspk_vec, time_mat

def load_data2(filebase, shank_num, name):
   #f = open(filebase+'/info')
   #info = f.read()
   #lst_info = info.split('\n')
   #session_name = lst_info[0]

    clu = np.load(filebase + '/' + name + '.clu.' + shank_num + '.npy')
    res = np.load(filebase + '/' + name + '.res.' + shank_num + '.npy')
    cc = np.load(filebase + '/' + name + '.cc.' + shank_num + '.npy')
    time_mat = np.load(filebase + '/' + name + '.timeMat.' + shank_num + '.npy')
    mspk = np.load(filebase + '/' + name + '.mspk.' + shank_num + '.npy')
    sspk = np.load(filebase + '/' + name + '.sspk.' + shank_num + '.npy')
    nspk_vec = np.load(filebase + '/' + name + '.nspk_vec.' + shank_num + '.npy')
    nspk_vec = np.squeeze(nspk_vec)
    return clu, res, cc, mspk, sspk, nspk_vec, time_mat


def clu_organize2(filebase, sorted_clu, idx_ignored, n_original,shank_num, name):
    new_clu = np.ones((n_original, 1))
    new_clu[idx_ignored] = 0
    c = 0
    for i in range(n_original):
        if new_clu[i][0] == 1:
            new_clu[i][0] = sorted_clu[c]
            c += 1
    fname = filebase + "/" + name + ".clu." + shank_num + "_2"
    np.save(fname, new_clu)
    return new_clu

def main_loop(clu, res, cc, mean_spk, std_spk, nspk_vec, time_mat):
    n_original = len(clu)

    idx_ignored = np.where((clu == 0) | (clu == 1))
    clu = np.delete(clu, idx_ignored[0], 0)

    #thresholds = [0.999, 0.995, 0.98, 0.7, 0.6, 0.5]
    #thresholds = [0.999, 0.8, 0.6]
    thresholds = [0.999, 0.995, 0.95, 0.9]
    unique_clu = np.unique(clu)
    for k in range(2):
        for i in thresholds:
            j = 0
            while j < len(unique_clu)-1:
                con = True
                while con:
                    feat_mat = np.empty((1, 5, 128))
                    n = len(unique_clu)

                    for m in range(j+1, n):
                        x = prepare_features_MC(j, m, clu, mean_spk, std_spk, cc, time_mat, unique_clu,k)
                        feat_mat = np.concatenate((feat_mat, x), 0)
                    if len(feat_mat) > 1:
                        feat_mat = np.delete(feat_mat, 0, axis=0)
                        a = predict2(feat_mat)
                        idx = np.where(a[:, 1] > i)
                        idx = idx[0]
                    else:
                        idx = np.empty(0)
                    if len(idx) < 1:
                        con = False
                        j += 1
                    else:
                        idx += (j+1)
                        if np.max(idx) <= len(unique_clu)-1:
                            g1 = unique_clu[idx]
                            g2 = np.array([unique_clu[j]])
                            group_units = np.concatenate((g1, g2))
                            clu, mean_spk, std_spk, nspk_vec, cc, time_mat = update_all2(cc, group_units, clu, mean_spk, std_spk, nspk_vec, time_mat, [0, 1])

                        unique_clu = np.unique(clu)

    new_clu = clu_organize(clu, idx_ignored, n_original)
    return new_clu


def clu_organize(sorted_clu, idx_ignored, n_original):
    new_clu = np.ones((n_original, 1))
    new_clu[idx_ignored] = 0
    c = 0
    for i in range(n_original):
        if new_clu[i][0] == 1:
            new_clu[i][0] = sorted_clu[c]
            c += 1
    # save_clu(filebase, new_clu, shank_num)
    return new_clu


def save_clu(filebase, new_clu, shank_num):
    f = open(filebase + '/info')
    info = f.read()
    lst_info = info.split('\n')
    session_name = lst_info[0]

    fname = filebase + "/" + session_name + ".clu." + shank_num + "_2"
    np.save(fname, new_clu)
    return


def getX(file):
    data = np.load(file)
    x = data[:, :128]
    x[4,0:83] += -0.5
    return torch.FloatTensor(x)


def getY(file):
    data = np.load(file)
    y = data[0, -1]
    y = str(int(y))
    return y

