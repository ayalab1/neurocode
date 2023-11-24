from fastai.vision.all import *
from inception_time import *
from fastai.callback import *
from fastai.text.all import *
import numpy as np
import torch
from prepare_features_nc import *
import os

def getX(item):
    x = X1[:, :, item].float()
    x[2, 0:41] += -0.5
    return x


def getY(item):
    y = str(int(Y[item]))
    return y


def predict(fet_mat):
    import pathlib
    posix_backup = pathlib.PosixPath
    try:
        pathlib.PosixPath = pathlib.WindowsPath
        learn = load_learner(os.path.join(pathlib.Path(__file__).parent.resolve(),"new_nc_301022.pkl"))
    finally:
        pathlib.PosixPath = posix_backup    

    dl = learn.dls.test_dl(fet_mat)
    preds = learn.get_preds(dl=dl)
    p = preds[0].numpy()
    p1 = np.argmax(p, axis=1)
    p1 = p1 + 1
    return p

def get_preds(clu, mean_spk, std_spk, cc, time_mat, u_clu):
    uclu = np.unique(clu)
    fet_mat = np.zeros((len(uclu),3,128))
    for i in range(len(uclu)):
        fet1 = prepare_features_nc(i, clu, mean_spk, std_spk, cc, time_mat, uclu)
        fet_mat[i,:,:]= fet1
    fet_mat =torch.tensor(fet_mat).float()
    preds = predict(fet_mat)
    return preds
    



