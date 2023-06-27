#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:21:12 2022

@author: lidor
"""
import numpy as np
from phylib.io.array import _clip
from phylib.utils import Bunch
from phylib.stats import ccg as CCG
import pandas as pd

def upsampCCG(lags2,lags,jit,cc1):
    upsampcc          =np.zeros_like(lags2) 
    for i in range(len(lags2)):
        T             = lags2[i]
        idx           = np.where([(lags>=T-jit) & (lags<T+jit)] )
        upsampcc[i]   = np.sum(cc1[idx[1]])
    return upsampcc

# compute CCG with 1 ms bin size
def compCCG(res,clu,FS=20000,window_size=0.05):
    res             = np.int64(np.squeeze(res))
    clu             = np.int64(np.squeeze(clu))
    
        
    Uclu            = np.unique(clu)
    BinSize         = 0.001
    ms              = 1000
    jit             = ms*BinSize/2
    
    if np.mod(window_size*ms, 2)!=0:
        window_size     = np.round(window_size*ms+1)/ms
    
    cc =CCG.correlograms(res/FS,clu,Uclu,sample_rate=FS,bin_size=0.0001,window_size=window_size)
    
    lags  = np.arange((-window_size/2)*ms,(window_size/2)*ms+0.0001*ms,0.1)
    lags2 = np.arange((-window_size/2)*ms,(window_size/2)*ms+BinSize*ms,BinSize*ms)
    Ncc   = np.zeros((len(Uclu),len(Uclu),len(lags2,)))
    
    for i in range(len(Uclu)):
        for j in range(len(Uclu)):
            cc1            = cc[i,j,:]
            Ncc[i,j,:]     = upsampCCG(lags2,lags,jit,cc1)
    Ncc            = np.transpose(Ncc, (2, 1, 0))
    return Ncc,lags2


# clu =np.load('mP79_16.clu.1.npy')
# res =np.load('mP79_16.res.1.npy')
# Ncc =compCCG(res,clu,20000,window_size=0.05)

# cluFname   = '/home/lidor/data/DBC/Code/AUSS_py/mP31_04.clu.2'
# resFname   = '/home/lidor/data/DBC/Code/AUSS_py/mP31_04.res.2'
#spkFname   = '/home/lidor/data/DBC/Code/AUSS_py/mP31_04.spk.2'
# clu        = pd.read_table(cluFname, delimiter=" ",header=None).to_numpy()
# clu        = np.delete(clu,0)
# res        = pd.read_table(resFname,delimiter=" ",header=None).to_numpy()
# res        = np.squeeze(res)
# Ncc        = compCCG(res,clu,FS=20000,window_size=0.040)



































   