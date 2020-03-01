# -*- coding: utf-8 -*-
"""
Created on Sun Mar 1 15:36:08 2020

@author: XuemeiGu
email: njuxmgu@gmail.com

"""

##### import packages 
import numpy as np
import sympy as sp
import random
from numpy.random import choice
from sympy.matrices import SparseMatrix  
from sympy import collect, expand, Symbol,sqrt,pi,I
from itertools import combinations, chain

FF1, FF2, FF3, FF4,FF5, HH, GG1, GG2, GG3, GG4,GG5=map(sp.IndexedBase,['FF1','FF2','FF3','FF4','FF5','HH','GG1','GG2','GG3','GG4','GG5'])  
l,l1, l2, l3, l4, l5, l6, l7, l8,coeff=map(sp.Wild,['l','l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8','coeff'])

zero=Symbol('zero') ## using in 4-fold coincidence 
psi=Symbol('psi')  ## represent a quantum state

## corresponding symbolic number 
sqr2=sqrt(2)/2  ## equals to 1/sqrt(2)
imagI=I ## equals to Imaginary
Pi=pi

  

## (*For calculating the Schmidt-Rank Vector*)
def toHH(expr):  ## change the state form to calculation the SRV
    expr1=expr.replace(coeff*FF2[l1]*FF3[l2]*FF4[l3]*FF5[l4],sp.conjugate(coeff)*GG2[l1]*GG3[l2]*GG4[l3]*GG5[l4])
    rho0=sp.expand(expr*expr1)
    rho0=rho0.replace(coeff*FF2[l1]*FF3[l2]*FF4[l3]*FF5[l4]*GG2[l5]*GG3[l6]*GG4[l7]*GG5[l8],coeff*HH[l1,l2,l3,l4,l5,l6,l7,l8])
    return rho0
 
def PartialTraceOne(expr,n): ##calculate the partial trace such as A|BCD 
    
    dictadd=collect(expr, [HH[l1,l2,l3,l4,l5,l6,l7,l8]], evaluate=False)
    TermsCoeff=list(dictadd.items())
    
    ParticleOne=[]
    ParticleTwo=[]
    ParticleThree=[]
    ## get the size of the matrix
    for ii in range(len(TermsCoeff)):
        HHList=TermsCoeff[ii][0]
        if HHList.indices[n-1]==HHList.indices[n+3]:
           ll=[HHList.indices[0],HHList.indices[1],HHList.indices[2],HHList.indices[3],HHList.indices[4],HHList.indices[5],HHList.indices[6],HHList.indices[7]]
           del(ll[n-1],ll[n+2])  ## because cannot del all at the same time, thus do it one by one, the index is not n+2
           
           ParticleOne.append(ll[0])
           ParticleTwo.append(ll[1])
           ParticleThree.append(ll[2])
    # start from 0      
    Upperone=max(ParticleOne)+1
    Lowerone=min(min(ParticleOne),0)
    Uppertwo=max(ParticleTwo)+1
    Lowertwo=min(min(ParticleTwo),0)
    Upperthree=max(ParticleThree)+1
    Lowerthree=min(min(ParticleThree),0)
    
    rangeP1=Upperone-Lowerone  
    rangeP2=Uppertwo-Lowertwo  
    rangeP3=Upperthree-Lowerthree 
    
    Msize=(rangeP1*rangeP2*rangeP3)
    SMatrix=SparseMatrix(Msize, Msize, {(0, 0): 0})
    
    for ii in range(len(TermsCoeff)):
        HHList=TermsCoeff[ii][0]
        if HHList.indices[n-1]==HHList.indices[n+3]:
           ll=[HHList.indices[0],HHList.indices[1],HHList.indices[2],HHList.indices[3],HHList.indices[4],HHList.indices[5],HHList.indices[6],HHList.indices[7]]
           del(ll[n-1],ll[n+2]) ## because cannot del all at the same time, thus do it one by one, the index is not n+2
           Dimrow=(ll[0]-Lowerone)*rangeP3*rangeP2+(ll[1]-Lowertwo)*rangeP3+(ll[2]-Lowerthree)
           Dimcol=(ll[3]-Lowerone)*rangeP3*rangeP2+(ll[4]-Lowertwo)*rangeP3+(ll[5]-Lowerthree)
           SMatrix=SparseMatrix(Msize, Msize, {(Dimrow,Dimcol):TermsCoeff[ii][1]})+SMatrix
    return SMatrix.rank() 


def PartialTraceTwo(expr,m,n): ##calculate the partial trace such as AB|CD 
    
    dictadd=collect(expr, [HH[l1,l2,l3,l4,l5,l6,l7,l8]], evaluate=False)
    TermsCoeff=list(dictadd.items())
    
    ParticleOne=[]
    ParticleTwo=[]
    ## get the size of the matrix
    for ii in range(len(TermsCoeff)):
        HHList=TermsCoeff[ii][0]
        if HHList.indices[m-1]==HHList.indices[m+3] and HHList.indices[n-1]==HHList.indices[n+3]:
           ll=[HHList.indices[0],HHList.indices[1],HHList.indices[2],HHList.indices[3],HHList.indices[4],HHList.indices[5],HHList.indices[6],HHList.indices[7]]
           del(ll[m-1]) 
           del(ll[m+2])  ## because cannot del all at the same time, thus do it one by one, the index is not n+2
           del(ll[n-2]) 
           del(ll[n]) 
           ParticleOne.append(ll[0])
           ParticleTwo.append(ll[1])
    # start from 0      
    Upperone=max(ParticleOne)+1
    Lowerone=min(min(ParticleOne),0)
    Uppertwo=max(ParticleTwo)+1
    Lowertwo=min(min(ParticleTwo),0)
    rangeP1=Upperone-Lowerone  
    rangeP2=Uppertwo-Lowertwo   
    
    Msize=(rangeP1*rangeP2)
    SMatrix=SparseMatrix(Msize, Msize, {(0, 0): 0})
    
    for ii in range(len(TermsCoeff)):
        HHList=TermsCoeff[ii][0]
        if HHList.indices[m-1]==HHList.indices[m+3] and HHList.indices[n-1]==HHList.indices[n+3] :
           ll=[HHList.indices[0],HHList.indices[1],HHList.indices[2],HHList.indices[3],HHList.indices[4],HHList.indices[5],HHList.indices[6],HHList.indices[7]]
          ## because cannot del all at the same time, thus do it one by one, the index is not n+2
           del(ll[m-1]) 
           del(ll[m+2])  ## because cannot del all at the same time, thus do it one by one, the index is not n+2
           del(ll[n-2]) 
           del(ll[n]) 
           Dimrow=(ll[0]-Lowerone)*rangeP2+(ll[1]-Lowertwo)
           Dimcol=(ll[2]-Lowerone)*rangeP2+(ll[3]-Lowertwo)
           SMatrix=SparseMatrix(Msize, Msize, {(Dimrow,Dimcol):TermsCoeff[ii][1]})+SMatrix
    return SMatrix.rank()

 
## The function to calculation the SRV of the quantum state
def SchmidtRankVector(expr):  
    rho=toHH(expr)
    if rho==0 :return [0,0,0,0,0,0,0]
    else: return [PartialTraceOne(rho,1),PartialTraceOne(rho,2),PartialTraceOne(rho,3),PartialTraceOne(rho,4),PartialTraceTwo(rho,1,2),PartialTraceTwo(rho,1,3),PartialTraceTwo(rho,1,4)]


## test: give a state (only consider same coefficient ) and returen the SRV number
States=FF2[0]*FF3[0]*FF4[0]*FF5[0]+FF2[1]*FF3[1]*FF4[1]*FF5[0]+FF2[0]*FF3[2]*FF4[3]*FF5[1];
SRV=SchmidtRankVector(States) 
print('SRV: ',SRV)
                        
