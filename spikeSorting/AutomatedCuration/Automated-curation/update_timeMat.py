import numpy as np

def update_time_mat(time_mat, group_units, clu, ignoredClu):
    unique_clu = np.unique(clu)
    unique_clu = np.setdiff1d(unique_clu, ignoredClu)
    n = len(group_units)
    x = np.where(unique_clu == group_units[0])
    new_mat = time_mat
    for i in range(n, 1, -1):
        y = np.where(unique_clu == group_units[i - 1])
        new_mat[x, :] = new_mat[x, :] + new_mat[y, :]
        new_mat = np.delete(new_mat, y, axis=0)
    return new_mat

# time_mat = np.load('/home/tali/data/mA234_35_all/npy_files/mA234_35_Shirly.timeMat.1.npy')
# clu = np.load('/home/tali/data/mA234_35_all/npy_files/mA234_35_Shirly.clu.1.npy')
# group_units = [7,15,40]
# update_time_mat(time_mat, group_units, clu, [0,1])