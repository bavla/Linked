
wdir = "C:/Users/batagelj/work/projects/BM/py"
gdir = "c:/users/batagelj/work/python/graph/Nets"
ddir = "C:/Users/batagelj/work/projects/BM/data"

import sys, os, re, json
sys.path = [gdir]+sys.path; os.chdir(wdir)
from Nets import Network as N
import numpy as np
from copy import copy, deepcopy
from timeit import default_timer as timer
from datetime import datetime

def savePajekPar(L,file):
    clu = open(file,'w');  n=len(L)
    clu.write('*vertices '+str(n)+'\n')
    for i in range(n): clu.write(str(L[i])+'\n')
    clu.close()

def normPar(L):
    I = -np.ones(1+max(L),dtype=int)
    N = np.zeros(len(L),dtype=int); k = -1
    for i,j in enumerate(L):
        if I[j]<0: k += 1; I[j] = k
        N[i] = I[j]
    return N

def P(net,C,Types):
    n = net._info['nNodes']; nC = Types.shape[0]; size = np.zeros(nC,dtype=int)
    for i in C: size[i] += 1
    dCnt = np.zeros(nC); Cnt = np.zeros((nC,nC))
  # Sum = np.zeros((nC,nC)); rSum = np.zeros((nC,n)); cSum = np.zeros((nC,n))
    rCnt = np.zeros((nC,n)); cCnt = np.zeros((nC,n))
    rPos = np.zeros((nC,nC)); cPos = np.zeros((nC,nC))
    error = np.zeros((nC,nC)); btype = np.zeros((nC,nC),dtype=int)
    for a in net.links():
        u,v,d,r,W = net.link(a); w = W['w']; p = u-1; q = v-1
        Cnt[C[p],C[q]] += 1; rCnt[C[q],p] += 1; cCnt[C[p],q] += 1
        if not d: Cnt[C[q],C[p]] += 1; rCnt[C[p],q] += 1; cCnt[C[q],p] += 1
        if p==q: dCnt[C[p]] += 1
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

def optimize(net,C,Types,minSize):
    nC = Types.shape[0]; size = np.zeros(nC,dtype=int)
    for i in C: size[i] += 1
    if np.min(size)<minSize: return(C,np.inf,None,None)
    found = True; Clu = np.array(C,dtype=int); Cm = Clu.copy()
    Err, Q, T = P(net,Clu,Types)
    while found:
        found = False
        for a in net.links():
            u,v,d,r,W = net.link(a); w = W['w']
            p = u-1; q = v-1; Cp = Clu[p]; Cq = Clu[q]
            if not(Cp==Cq):
                Clu[p] = Cq; size[Cp] += -1; size[Cq] += 1
                if size[Cp]<minSize: err = np.inf
                else: err, Q, T = P(net,Clu,Types)
                if err<Err: Cm = Clu.copy(); Err = err; found = True
                else:
                    Clu[p] = Cp; Clu[q] = Cp; size[Cp] += 2; size[Cq] += -2
                    if size[Cq]<minSize: err = np.inf
                    else: err, Q, T = P(net,Clu,Types)
                    if err<Err: Cm = Clu.copy(); Err = err; found = True
                    else:
                        Clu[p] = Cq; size[Cp] += -1; size[Cq] += 1
                        err, Q, T = P(net,Clu,Types)
                        if err<Err: Cm = Clu.copy(); Err = err; found = True
                        else: Clu[p] = Cp; Clu[q] = Cq 
    return (Cm, Err)

def searchBM(Types,nSteps=10,Pbest=np.inf,R=[1,1,1],Q=None,minSize=1,
             normalize=True,trace=1):
    print('netsBM:', datetime.now())
    start = timer(); n = net._info['nNodes']; nC = Types.shape[0]
    iD = np.array([net.inDegree(u) for u in net.nodes()],dtype='I')        
    oD = np.array([net.outDegree(u) for u in net.nodes()],dtype='I')
    ind = np.logical_and(iD<2,oD<2); nSpec = np.sum(ind); spec = np.any(ind)
    if Q is None: Q = [1]*nC
    r = np.array(R)/sum(R); q = np.array(Q)/sum(Q)
    found = False; bad = 0
    for k in range(nSteps):
        j = np.random.choice([0,1,2],p=r)
        if j==0: C = np.random.choice(range(nC),size=(n))
        elif j==1: C = np.random.choice(range(nC),p=q,size=(n))
        elif j==2:
            C = spec+np.random.choice(range(nC-spec),size=(n))
            if spec: C[ind] = 0
        R = optimize(net,C,Types,minSize); Er = R[1]
        if trace>1: print(k+1,j,Er)
        if Er == np.inf: bad += 1
        elif Er <= Pbest:
            if Er < Pbest: Pbest = Er; Cbest = []; nBest = 0
            found = True; nBest += 1; newC = True
            CC = normPar(R[0]) if normalize else R[0]
            for CB in Cbest:
                if np.all(CC==CB): newC = False; break
            if newC: Cbest.append(CC)
            if trace>0:
                if n > 50: print(k+1,j,Er)
                else: print(k+1,j,Er,"\n  C :",C,"\n  C*:",CC)
    end = timer(); print("time:",round(end-start,2),"s")
    if bad>0: print(bad,"bad initial partitions")
    if found: return (Pbest,nBest,Cbest)

# net = N.loadPajek(ddir+"/kite.net"); net.Info()
net = N.loadPajek(ddir+"/class.net"); net.Info()
# net = N.loadPajek(ddir+"/dolphins.net"); net.Info()
# net = N.loadPajek(ddir+"/USAir97.net"); net.Info()

XXX = 0; NUL = 1; COM = 2; REG = 3; RRE = 4; CRE = 5
n = net._info['nNodes']; nC = 3; Q = [3,3,9] 
strEq = np.array([{NUL,COM}]*(nC*nC),dtype='O').reshape(nC,nC)

# n = net._info['nNodes']; nC = 6
# Rez = searchBM(strEq,100,Q=Q)
# Pbest = Rez[0]; nBest = Rez[1]
# regEq = np.array([{NUL,COM,REG}]*(nC*nC),dtype='O').reshape(nC,nC)
# Reg = searchBM(regEq,30,Q=Q,minSize=2)
# Prbest = Reg[0]; nrBest = Reg[1]
# savePajekPar(Rez[2][0],"USairStr6.clu")
