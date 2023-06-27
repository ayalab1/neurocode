import numpy as np


def update_cc(cc, group_units, clu=[]):

    if clu == []:
        all_units = np.array(group_units)
    else:
        all_units = np.unique(clu)

    x = np.where(all_units == group_units[0])
    new_cc = cc
    n = len(group_units)
    for i in range(n-1, 0, -1):
        idx = np.where(all_units == group_units[i])
        new_cc[:, x, :] = new_cc[:, x, :] + new_cc[:, idx, :]
        new_cc = np.delete(new_cc, idx, axis=1)

    for i in range(n-1, 0, -1):
        idx = np.where(all_units == group_units[i])
        new_cc[:, :, x] = new_cc[:, :, x] + new_cc[:, :, idx]
        new_cc = np.delete(new_cc, idx, axis=2)

    return new_cc

#clu = np.load('/home/tali/matlab/AUSS_python/mP31_04.clu.1.npy')
#cc = np.load('/home/tali/matlab/AUSS_python/mP31_04.cc.1.npy')
#group_units = [6, 15]

#new = update_cc(cc, group_units)


