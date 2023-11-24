from update_cc import *
from update_mspk import *
from update_timeMat import *


def update_all(cc, group_units, clu, mspk, sspk, nspk_vec, ignored_clu=[]):
    new_cc = cc
    new_clu = clu
    new_mspk = mspk
    new_sspk = sspk
    new_spk_vec = nspk_vec

# Removes 0 from clu, maybe move to another function
#    idx = np.where(new_clu == 0)
#    new_clu = np.delete(new_clu, idx)
    group_units.sort()
    new_mspk, new_sspk, new_spk_vec = update_mspk(new_mspk, new_sspk, new_spk_vec, group_units, new_clu, ignored_clu)
    new_cc = update_cc(new_cc, group_units, new_clu)
    new_clu = update_clu(new_clu, group_units)

    return new_clu, new_mspk, new_sspk, new_spk_vec, new_cc


def update_clu(clu, group_units):
    new_clu = clu
    n = len(group_units)
    for j in range(1, n):
        idx = np.where(new_clu == group_units[j])
        new_clu[idx] = group_units[0]
    return new_clu


def update_all2(cc, group_units, clu, mspk, sspk, nspk_vec, time_mat, ignored_clu=[]):
    group_units.sort()
    mspk, sspk, nspk_vec = update_mspk(mspk, sspk, nspk_vec, group_units, clu, ignored_clu)
    cc = update_cc(cc, group_units, clu)
    time_mat = update_time_mat(time_mat, group_units, clu, ignored_clu)
    clu = update_clu(clu, group_units)
    return clu, mspk, sspk, nspk_vec, cc, time_mat


#clu = np.load('/home/tali/matlab/AUSS_python/mP31_04.clu.1.npy')
#cc = np.load('/home/tali/matlab/AUSS_python/mP31_04.cc.1.npy')
#mspk = np.load('/home/tali/matlab/AUSS_python/mP31_04.mspk.1.npy')
#sspk = np.load('/home/tali/matlab/AUSS_python/mP31_04.sspk.1.npy')
#nspk_vec = np.load('/home/tali/matlab/AUSS_python/mP31_04.nspk_vec.1.npy')
#nspk_vec = np.squeeze(nspk_vec)
#group_units = [3, 4, 5, 29, 6]
#update_all(cc, group_units, clu, mspk, sspk, nspk_vec, [0,1])

#Nmspk = np.load('/home/tali/matlab/AUSS_python/mP31_04/mP31_04.Nmspk_t.npy')
#Nsspk = np.load('/home/tali/matlab/AUSS_python/mP31_04/mP31_04.Nsspk_t.npy')
#Nspk_vec = np.load('/home/tali/matlab/AUSS_python/mP31_04/mP31_04.Nspk_vec_t.npy')