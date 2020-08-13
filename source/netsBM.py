
wdir = "C:/Users/batagelj/work/projects/BM/py"
gdir = "c:/users/batagelj/work/python/graph/Nets"
ddir = "C:/Users/batagelj/work/projects/BM/data"

import sys, os, re, datetime, json
sys.path = [gdir]+sys.path; os.chdir(wdir)
from TQ import *
from Nets import Network as N
import numpy as np
from copy import copy, deepcopy
# from hcluRCpq import printNet

def P(net,C,Types):
    n = net._info['nNodes']; nC = 1+np.max(C); size = np.zeros(nC)
    for i in C: size[i] += 1
    dCnt = np.zeros(nC); Cnt = np.zeros((nC,nC))
  # Sum = np.zeros((nC,nC))
  # rSum = np.zeros((nC,n)); cSum = np.zeros((nC,n))
    rCnt = np.zeros((nC,n)); cCnt = np.zeros((nC,n))
    rPos = np.zeros((nC,nC)); cPos = np.zeros((nC,nC))
    error = np.zeros((nC,nC)); btype = np.zeros((nC,nC),dtype=int)
    for a in net.links():
        u,v,d,r,W = net.link(a); w = W['w']; p = u-1; q = v-1
        Cnt[C[p],C[q]] += 1
        rCnt[C[q],p] += 1; cCnt[C[p],q] += 1
        if p==q:
            dCnt[C[p]] += 1
    rB = (rCnt>0).astype(int); cB = (cCnt>0).astype(int)        
    for i in range(nC):
        for s in range(n):
            rPos[i,C[s]] += rB[i,s]; cPos[C[s],i] += cB[i,s]
    for i in range(nC):
        for j in range(nC):
            err = np.inf; typ = XXX
            if NUL in Types[i,j]: 
                ert = Cnt[i,j]
                if i==j: ert += min(0,size[i]-2*dCnt[i])
                if ert < err: err = ert; typ = NUL
            if COM in Types[i,j]: 
                ert = size[i]*size[j] - Cnt[i,j]
                if i==j: ert += min(0,2*dCnt[i]-size[i])
                if ert < err: err = ert; typ = COM
            if REG in Types[i,j]: 
                ert = (size[j]-cPos[i,j])*size[i] +\
                      (size[i]-rPos[i,j])*cPos[i,j]
                if ert < err: err = ert; typ = REG
            if RRE in Types[i,j]: 
                ert = (size[i]-rPos[i,j])*size[j]
                if ert < err: err = ert; typ = RRE
            if CRE in Types[i,j]: 
                ert = (size[j]-cPos[i,j])*size[i] 
                if ert < err: err = ert; typ = CRE
            error[i,j] = err; btype[i,j] = typ
    return (np.sum(error), error, btype)

def optimize(net,C,Types):
    found = True; Clu = np.array(C,dtype=int); Cm = Clu.copy()
    Err, Qm, Tm = P(net,Clu,Types)
    while found:
        found = False
        for a in net.arcs():
            u,v,d,r,W = net.link(a); w = W['w']
            p = u-1; q = v-1; Cp = Clu[p]; Cq = Clu[q]
            if not(Cp==Cq):
                Clu[p] = Cq
                err, Q, T = P(net,Clu,Types)
                if err<Err:
                    Cm = Clu.copy(); Err = err; found = True
                    Qm = Q.copy(); Tm = T.copy()
                else:
                    Clu[p] = Cp; Clu[q] = Cp
                    err, Q, T = P(net,Clu,Types)
                    if err<Err:
                        Cm = Clu.copy(); Err = err; found = True
                        Qm = Q.copy(); Tm = T.copy()
                    else:
                        Clu[p] = Cq
                        err, Q, T = P(net,Clu,Types)
                        if err<Err:
                            Cm = Clu.copy(); Err = err; found = True
                            Qm = Q.copy(); Tm = T.copy()
    return (Cm, Err, Qm, Tm)

def searchBM(nSteps=10,Pbest=np.inf):
    found = False
    for k in range(nSteps):
        j = k % 3
        if j==0: C = (nC*rng.random(n)).astype(int)
        elif j==1: C = C0.copy(); rng.shuffle(C)
        elif j==2:
            iD = np.array([net.inDegree(u) for u in net.nodes()],dtype='I')        
            oD = np.array([net.outDegree(u) for u in net.nodes()],dtype='I')
            C = 1+((nC-1)*rng.random(n)).astype(int)
            C[np.logical_and(iD<2,oD<2)] = 0
        R = optimize(net,C,Types); Er = R[1]
        if Er < Pbest:
            Pbest = Er; Rbest = deepcopy(R); found = True
            print(k+1,j,Er," C*:",R[0])
    if found: return Rbest

# net = N.loadPajek(ddir+"/kite.net"); net.Info()

ddir = "C:/Users/batagelj/work/Python/graph/Nets/BM/dat"
net = N.loadPajek(ddir+"/class.net"); net.Info()

XXX = 0; NUL = 1; COM = 2; REG = 3; RRE = 4; CRE = 5
rng = np.random.default_rng()
n = net._info['nNodes']; nC = 3
Types = np.array([{NUL,COM}]*(nC*nC),dtype='O').reshape(nC,nC)
C0 = np.array([0]*3 + [1]*3 + [2]*9, dtype='I')
# C0 = np.array([ i%nC for i in range(n) ], dtype='I')
# Rez = searchBM(100)
# Pbest = Rez[1]; Rbest = deepcopy(R)
# Rez = searchBM(100,Pbest)
