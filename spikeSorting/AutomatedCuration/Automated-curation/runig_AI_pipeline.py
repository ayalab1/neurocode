#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 15:12:11 2022

@author: lidor
"""
import numpy as np
from ccg import *
import create_feat_tsc as tsc
import sort_shank as MC
from loading_data import *
from noise_classifier import get_preds

def getX(item):
    x = X1[:,item].float().unsqueeze(0)
    return x

def getY(item):
    y = str(int(Y[item]))
    return y

# inputs 
filebase   = '/home/lidor/data/AUSS_project/Automated_curation/test_data/mP31_04'
shank      = 2
Nchannels  = 10
Nsamp      = 32

# Loading data
clu,res    = load_clures(filebase,shank)
mspk,sspk  = make_mSPK(filebase,shank,Nchannels=Nchannels,Nsamp=Nsamp)
cc         = get_CCmat(filebase,shank)
u_clu      = np.unique(clu)

# generating time mat for NC
time_mat   = get_time_mat1(res, clu)

# get porpabilities form NC
pred     = get_preds(clu, mspk, sspk, cc, time_mat, u_clu)

# generat a clu which noise and multuybuts are lebeled as zero
cleanClu = tsc.tsc(pred,clu)
ind,Z    = tsc.get_cleanClu_idx(clu,cleanClu)

# make new featurs for the MC
nspk_vec           = tsc.compute_Nvec(cleanClu)[1:]
cluster_ids        = np.unique(cleanClu)
time_mat           = tsc.compute_timeMat(cleanClu,res,cluster_ids)[1:,:]
mean_spk,std_spk   = tsc.orgnize_WF(mspk,sspk,ind,Z)

sample_rate    = 20000
cc             = compCCG(res,cleanClu,FS=sample_rate,window_size=0.042)[0]
cc             = cc[1:-1,1:,1:]

newCLu         = MC.main_loop(cleanClu, res, cc, mean_spk, std_spk, nspk_vec, time_mat)

clu1        = clu
clu2        = newCLu
U_id        = np.unique(clu2)
reco_list   = list()
for i in U_id:
    idx = np.where(clu2==i)[0]
    l   = np.unique(clu1[idx])
    if len(l) > 1:
        reco_list.append(l)
print(reco_list)
