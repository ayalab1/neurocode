#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 12:30:57 2022

@author: lidor
"""
import numpy as np
import pandas as pd
from ccg import *


def load_clures(filebase,shank):
    shank     = str(shank)
    clu_ex    = '.clu.'
    res_ex    = '.res.'
    cluFname  =filebase+clu_ex+shank 
    resFname  = filebase+res_ex+shank

    clu        = pd.read_table(cluFname, delimiter=" ",header=None).to_numpy()
    clu        = np.delete(clu,0)
    res        = pd.read_table(resFname, delimiter=" ",header=None).to_numpy()
    res        = np.squeeze(res)
    return clu,res


def make_mSPK(filebase,shank,Nchannels=10,Nsamp= 32):
    clu,res   = load_clures(filebase,shank)
    shank     = str(shank)
    spk_ex    = '.spk.'
    spkFname  = filebase+spk_ex+shank
    
    
    
    Nvar               = Nsamp*Nchannels
    Nspkblock          = 10000
    Nspk               = len(res)
    vec                = np.append(np.arange(0,Nspk,Nspkblock),Nspk)
    vec                = vec.astype(np.int64)
    Uclu               = np.unique(clu)
    mSPK               = np.zeros((len(Uclu),Nsamp,Nchannels))
    sSPK               = np.zeros((len(Uclu),Nsamp,Nchannels))
    nSPK               =np.zeros_like(Uclu) 
    for i,j in enumerate(vec[:-1]):
        dif             = vec[i+1]-j
        count           = dif*Nchannels*Nsamp
        ofset           = j*Nchannels*Nsamp
        ofset = ofset.astype(np.int64) # to prevent overflow in next line
        # memory map the file and select the relevant part according to the offset and count
        spk = np.memmap(spkFname, dtype=np.int16, mode="r", offset=ofset*2)[:count]
        spk             = np.reshape(spk,(dif,Nsamp,Nchannels),order='C')
        spk             = np.int64(spk)
        A               = int(ofset/Nvar)
        B               = int((ofset+count)/Nvar)
        clu1            = clu[A:B]
        uclu1           = np.unique(clu1)
        for zz in range(len(uclu1)):
            ci                = np.where(Uclu==uclu1[zz])[0]
            idx1              = np.where(clu1==uclu1[zz])[0]
            mSPK[ci,:,:]      = mSPK[ci,:,:] + np.expand_dims(np.sum(spk[idx1,:,:],axis=0),0) 
            sSPK[ci,:,:]      = sSPK[ci,:,:] + np.expand_dims(np.sum(spk[idx1,:,:]**2,axis=0),0)
            nSPK[ci]          = nSPK[ci]+np.sum(clu1==uclu1[zz])
        
    
    mspk  = mSPK/np.expand_dims(nSPK,(1,2))
    sspk  = np.sqrt( sSPK/np.expand_dims(nSPK,(1,2)) - mspk**2 )
    mspk   = np.transpose(mspk,(2,1,0))
    sspk   = np.transpose(sspk,(2,1,0))
    
    return(mspk,sspk)


def get_CCmat(filebase,shank):
    clu,res  = load_clures(filebase,shank)
    cc       = compCCG(res,clu,FS=20000,window_size=0.042)
    cc       = cc[0][1:-1,:,:]
    return cc

def get_time_mat1(res, clu):
    if len(res) == 0:
        T = []
    else:
        timeVec          = np.linspace(res[0],res[-1],num=88)
        timeVec          = np.int64(timeVec)
        uClu             = np.unique(clu)
        T                = np.zeros((len(uClu),87))
        for i in range(len(uClu)):
            uName        = uClu[i]
            idx          = clu == uName
            t1           = res[idx]
            v1           = np.zeros((87,1))
            for k in range(87):
                start    = timeVec[k]
                end1     = timeVec[k+1]
                n1       = len(t1[(t1>start) & (t1<end1)])
                v1[k]    = n1
            m1           = np.mean(v1)
            threshold1   = m1*0.1
            v1[v1<=threshold1] = 0
            v3           = np.zeros((87,1))

            for k in range(87):
                if v1[k] == 0:
                    v3[k] = -0.5
                else:
                    v3[k] = 0.5
            T[i,:]       = np.squeeze(v3)
    return T
    
    
if __name__ == '__main__':    
    #filebase   = r'D:\autocluster\V1JeanDidier_230619'
    #shank      = 2
    filebase = 'D:\\autocluster\\V1JeanDidier_230619\\V1JeanDidier_230619'
    shank = 3
    clu,res    = load_clures(filebase,shank)
    mspk,sspk  = make_mSPK(filebase,shank,Nchannels=16)
    cc         =get_CCmat(filebase,shank)
    time_mat   = get_time_mat1(res, clu)

