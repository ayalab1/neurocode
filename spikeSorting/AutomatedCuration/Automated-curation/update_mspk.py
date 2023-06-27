import numpy as np


def update_mspk(mspk, sspk, nspk_vec, group_units, clu, ignoredClu=[]):
    unique_clu = np.unique(clu)
    unique_clu = np.setdiff1d(unique_clu, ignoredClu)
    n = len(group_units)
    x = np.where(unique_clu == group_units[0])

    sum_spk = np.squeeze(mspk[:, :, x] * nspk_vec[x])
    new_mspk = mspk
    new_nspk_vec = nspk_vec
    new_sspk = sspk
    for i in range(n, 1, -1):
        mean1 = sum_spk / new_nspk_vec[x]
        n1 = new_nspk_vec[x]

        idx = np.where(unique_clu == group_units[i-1])

        nspk = new_nspk_vec[idx]
        sum_spk = sum_spk + (np.squeeze(new_mspk[:, :, idx] * nspk))
        new_nspk_vec[x] = new_nspk_vec[x] + new_nspk_vec[idx]

        n2 = new_nspk_vec[idx]
        mean2 = np.squeeze(mspk[:, :, idx])
        mean12 = sum_spk / new_nspk_vec[x]
        d1 = mean1 - mean12
        d2 = mean2 - mean12
        s1 = np.squeeze(new_sspk[:, :, x])
        s2 = np.squeeze(sspk[:, :, idx])

        tmp = np.expand_dims(np.sqrt((n1 * (s1 ** 2 + d1 ** 2) + n2 * (s2 ** 2 + d2 ** 2)) / (n1 + n2)), axis=2)
        new_sspk[:, :, x] = np.expand_dims(tmp, axis=3)

        new_nspk_vec = np.delete(new_nspk_vec, idx)
        new_mspk = np.delete(new_mspk, idx, axis=2)
        new_sspk = np.delete(new_sspk, idx, axis=2)

    tmp = np.expand_dims(sum_spk / new_nspk_vec[x], axis=2)
    new_mspk[:, :, x] = np.expand_dims(tmp, axis=3)

    return new_mspk, new_sspk, new_nspk_vec


#clu = np.load('/home/tali/matlab/AUSS_python/mP31_04.clu.1.npy')
#mspk = np.load('/home/tali/matlab/AUSS_python/mP31_04.mspk.1.npy')
#sspk = np.load('/home/tali/matlab/AUSS_python/mP31_04.sspk.1.npy')
#nspk_vec = np.load('/home/tali/matlab/AUSS_python/mP31_04.nspk_vec.1.npy')
#nspk_vec = np.squeeze(nspk_vec)
#group_units = [6, 15]

#f = update_mspk(mspk, sspk, nspk_vec, group_units, clu, [0, 1])