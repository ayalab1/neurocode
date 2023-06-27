#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:14:04 2022

@author: lidor
"""
import numpy as np

def tsc(pred,cluS):
 
    m        = np.max(pred,1)
    i        = np.argmax(pred,1)
    pred2    = pred
    
    Npred    = np.zeros((pred.shape[0],1))
    idx      = i==2
    Npred[idx,0]   = pred2[idx,2]+((pred2[idx,2]-pred2[idx,0]))+10
    idx          = i==1
    idx2         = np.logical_and(idx ,pred2[:,2]>=pred2[:,0])
    Npred[idx2,0]  = pred2[idx2,1]+5 +(pred2[idx2,2]-pred2[idx2,0]) 
    
    idx3           = np.logical_and(idx,  pred2[:,2]<pred2[:,0])
    Npred[idx3,0]  = pred2[idx3,1]+2 +(pred2[idx3,0]-pred2[idx3,2])
    
    idx          = i==0
    #Npred[idx,0]   =1-( pred2[idx,0]+(pred2[idx,0]-pred2[idx,2]))
    Npred[idx,0]   =1-( pred2[idx,0])
    pred         = Npred
    # orgenize new clu file
    sortPred     = np.sort(pred,axis=0)[::-1]
    idxP         = np.argsort(pred,axis=0)[::-1]
    newClu       = np.zeros_like(cluS)
    cleanClu     = np.zeros_like(cluS)
    u            = np.unique(cluS)
    nclu         = len(u)
    for i  in range(nclu):
         ind          = cluS == u[idxP[i,0]]
         newClu[ind]  = i+2
         if pred[idxP[i,0]]>0.01:
             cleanClu[ind]  = i+2
         elif pred[idxP[i,0]]<=0.01:
             cleanClu[ind]  = 0
    
    return cleanClu
         
def get_cleanClu_idx(clu,cleanClu):
    clu1    = np.int32(clu)
    clu2    = np.int32(cleanClu)
    Uclu1   = np.unique(clu1)
    Uclu2   = np.unique(clu2)
    ind     = np.zeros_like(Uclu1)
    Z       = np.zeros_like(Uclu1)>1
    C       = 0
    for i,u1 in enumerate(Uclu2):
        idx    = clu2==u1
        idx2   = np.isin(Uclu1,np.unique(clu1[idx]))
        l      = np.where(idx2)[0]
        ind[C:len(l)+C] = l
        Z[C:len(l)+C]   = np.unique(clu2[idx])>1
        C               = C+len(l)
    return ind,Z
      

def orgnize_WF(mspk,sspk,ind,Z):
    new_mspk    = mspk[:,:,ind]
    new_mspk    = new_mspk[:,:,Z]
    
    new_sspk    = sspk[:,:,ind]
    new_sspk    = new_sspk[:,:,Z]
    return new_mspk,new_sspk

def compute_timeMat(clu,res,cluster_ids):
    div_fac = 44
    timeVec = np.linspace(res[0],res[len(res)-1],div_fac+1,dtype="int32")
    T = np.zeros((len(cluster_ids),div_fac))
    for index, unit in enumerate(cluster_ids):
        idx = clu == unit
        t1    = res[idx] 
        v1   = np.zeros((div_fac))
        for k in range(div_fac):
            start    = timeVec[k]
            end1     = timeVec[k+1]
            n1       = len(  np.where(((t1>start) & (t1<end1) ) )[0] )
            v1[k]    = n1
        T[index]   = v1        
    return T

def compute_Nvec(clu):
    uni , counts = np.unique(clu, return_counts=True)
    return counts

































