import numpy as np


def prepare_features_nc(i, clu, mean_spk, std_spk, cc, time_mat, u_clu):
    mean_spk1 = mean_spk[:, :, i]
    [mean_spk1, ind] = trim_spk_4ch(mean_spk1)
    mean_spk1 = mean_spk1.flatten()
    max1 = max(abs(mean_spk1))

    std_spk1 = (std_spk[ind, :, i]).flatten()

    if np.sum(cc[:, i, i]) == 0:
        acc1 = cc[:, i, i]
    else:
        acc1 = cc[:, i, i] / np.max(cc[:, i, i])

    t = time_mat[i,:]

    last_row = np.concatenate((acc1.T, t.flatten()))
    last_row[0:41]+=-0.5
    x = (np.concatenate((mean_spk1, std_spk1))) / max1
    x = np.concatenate((x, last_row))
    x = np.reshape(x, (3, 128, 1))
    x = np.moveaxis(x, -1, 0)

    return x


def trim_spk_4ch(mean_spk):
    n_channels = np.size(mean_spk, 0)

    if n_channels < 4:
        new_mspk = mean_spk
        channels_idx = np.arange(0, n_channels)
        channels_idx = channels_idx.T
    else:
        M1 = np.amax(mean_spk, axis=1)
        M2 = np.amin(mean_spk, axis=1)
        I = np.argsort(np.abs(M1 - M2))
        channels_idx = I[-4:]
        channels_idx = np.flip(channels_idx)
        new_mspk = mean_spk[channels_idx.T, :]

    return new_mspk, channels_idx

