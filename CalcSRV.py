# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 23:36:08 2019

@author: XuemeiGu
email: njuxmgu@gmail.com

"""

##### import packages 
import numpy as np
import sympy as sp

from sympy.matrices import SparseMatrix  
from sympy import collect, expand, Symbol,sqrt,pi,I

 
a, b, c, d, e, f, FF1, FF2, FF3, FF4, HH, GG1, GG2, GG3, GG4=map(sp.IndexedBase,['a','b','c','d','e', 'f', 'FF1','FF2','FF3','FF4','HH','GG1','GG2','GG3','GG4'])  
l,l1, l2, l3, l4, l5, l6, l7, l8, x1, x2, x3, x4, x5, x6, coeff, powern =map(sp.Wild,['l','l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'coeff', 'powern'])

zero=Symbol('zero') ## using in 4-fold coincidence 
psi=Symbol('psi')  ## represent a quantum state

## corresponding symbolic number 
sqr2=sqrt(2)/2  ## equals to 1/sqrt(2)
imagI=I ## equals to Imaginary
Pi=pi

 
l_list=np.arange(-8,8,1).tolist() ## l_list define the dimension of OAM I used.
lnum=len(l_list)

## SPDC process
def DownConvOAM(lorder,p1,p2):  ## create the initial state,DC is a parameter
    initial_state=0
    for ii in range(-lorder,lorder+1):
        initial_state=p1[ii]*p2[-ii]+initial_state
    return initial_state

## define functions for OAM modes 
  
def BS_fun(expr,p1,p2): 
    if expr.base==p1:return expr.replace(p1[l],sqr2*(p2[l]+imagI*p1[-l]),map=False, simultaneous=True, exact=False)
    else: return expr.replace(p2[l],sqr2*(p1[l]+imagI*p2[-l]),map=False, simultaneous=True, exact=False)
    
def BS(psi,p1,p2):
    psi0=sp.expand(psi.replace(lambda expr: expr.base in [p1,p2], lambda expr: BS_fun(expr,p1,p2)))
    return psi0

def LI_fun(expr,p1,p2): 
    if expr.base==p1:return expr.replace(p1[l1],(sp.cos(l1*Pi/2)**2)*p1[l1]+imagI*(sp.sin(l1*Pi/2)**2)*p2[-l1],map=False, simultaneous=True, exact=False)
    else: return expr.replace(p2[l1],-(sp.cos(l1*Pi/2)**2)*p2[l1]+imagI*(sp.sin(l1*Pi/2)**2)*p1[-l1],map=False, simultaneous=True, exact=False)

def LI(psi,p1,p2):
    psi0=sp.expand(psi.replace(lambda expr: expr.base in [p1,p2], lambda expr: LI_fun(expr,p1,p2)))
    return psi0
  
def Reflection(expr, p):
    expr=imagI*expr.replace(p[l1],p[-l1], map=False, simultaneous=True, exact=False)
    return expr

def OAMHolo(expr, p, n):
    expr=expr.replace(p[l1],p[l1+n], map=False, simultaneous=True, exact=False)
    return expr

def DP(expr, p):
    expr=expr.replace(p[l1],imagI*sp.exp(imagI*l1*(Pi))*p[-l1], map=False, simultaneous=True, exact=False)
    return expr
 
## make 4-fold coincidence
def replaceRule(expr, repls):
    for k, m in repls.items():
        expr = expr.replace(k, m, map=False, simultaneous=True, exact=False)   
    return expr

def MakeFF2(expr): 
    NFoldrepls = {coeff*a[l1]*a[l2] : 0, coeff*b[l1]*b[l2] : 0, coeff*c[l1]*c[l2] : 0, coeff*d[l1]*d[l2] : 0}
    expr1=replaceRule(expr,NFoldrepls)
    NFoldrepls={coeff*a[l2]*b[l3]*c[l4]*d[l5]: 0}
    expr2=replaceRule(expr1,NFoldrepls)
    expr=expr1-expr2
    expr=expr.replace(coeff*a[l1]*b[l2]*c[l3]*d[l4],coeff*FF1[l1]*FF2[l2]*FF3[l3]*FF4[l4])
    return expr

## trigger
def Trigger(expr,p1,nlist): ##must be used after MakeFF2 in path p1
    lt=sp.Wild('lt',exclude=nlist)
    expr1=expr.replace(p1[lt],0)
    expr=expr1.replace(p1[l],1)
    return expr

## functions for maximum entangled SRV state
    
# retutn [(terms, coeff),(terms, coeff),...]    
def TermsCoeffList(expr): # retutn [(terms, coeff),(terms, coeff),...]
    
    dictadd=collect(expr, [FF2[x1]*FF3[x2]*FF4[x3]], evaluate=False)
    TermsCoeff=list(dictadd.items())
    return TermsCoeff

## return the number of the terms in the 4-fold state and check the equality for coeffs 
def CheckTermCoeffOfState(TermsCoeff): 
     
    coeffvalue={}
    for ii in range(len(TermsCoeff)):
        coefftemp=TermsCoeff[ii][1]
#        print('coefftemp   ',coefftemp)
        coeffvalue[ii]=expand(coefftemp*sp.conjugate(coefftemp))     
    Equlflag=0
    for iii in range(len(TermsCoeff)):
        if iii < len(TermsCoeff)-1:
            if coeffvalue[iii]==coeffvalue[iii+1]:
               Equlflag=1
            else:
               Equlflag=0
               break
    return [len(TermsCoeff), Equlflag]   


def AllTypesOfFFl(TermsCoeff,p1,n): # return the list for dimention in FFl (FF2 FF3 FF4)
    
    dictFFl={}
    dictFFn={}
    for ii in range(len(TermsCoeff)):
        dictFFl[ii]=TermsCoeff[ii][0]
        dictFFltemp=dictFFl[ii].as_ordered_factors() # decompose into something as [FF2[0], FF3[0], FF4[0]]
        dictFFn[ii]=dictFFltemp[n]
        
    FFldimtemp=list(dictFFn.values())  ## put all the corresponding term in a list
    TypesOfFF1=set(FFldimtemp)  ## remove repeated terms 
    return TypesOfFF1 # return the list for dimention in FFl (FF2 FF3 FF4)


## (*For calculating the Schmidt-Rank Vector*)
def toHH(expr):  ## change the state form to calculation the SRV
    expr1=expr.replace(coeff*FF2[l1]*FF3[l2]*FF4[l3],sp.conjugate(coeff)*GG2[l1]*GG3[l2]*GG4[l3])
    rho0=sp.expand(expr*expr1)
    rho0=rho0.replace(coeff*FF2[l1]*FF3[l2]*FF4[l3]*GG2[l4]*GG3[l5]*GG4[l6],coeff*HH[l1,l2,l3,l4,l5,l6])
    return rho0
 
def PartialTrace(expr,n): ##calculate the partial trace of the n_th photon 
    
    dictadd=collect(expr, [HH[l1,l2,l3,l4,l5,l6]], evaluate=False)
    TermsCoeff=list(dictadd.items())
    
    ParticleOne=[]
    ParticleTwo=[]
    ## get the size of the matrix
    for ii in range(len(TermsCoeff)):
        HHList=TermsCoeff[ii][0]
        if HHList.indices[n-1]==HHList.indices[n+2] :
           ll=[HHList.indices[0],HHList.indices[1],HHList.indices[2],HHList.indices[3],HHList.indices[4],HHList.indices[5]]
           del(ll[n-1],ll[n+1])  ## because cannot del all at the same time, thus do it one by one, the index is not n+2
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
        if HHList.indices[n-1]==HHList.indices[n+2] :
           ll=[HHList.indices[0],HHList.indices[1],HHList.indices[2],HHList.indices[3],HHList.indices[4],HHList.indices[5]]
           del(ll[n-1],ll[n+1])  ## because cannot del all at the same time, thus do it one by one, the index is not n+2
    #       print('rest: ',ll)
    #       print('rest: ',ll[0]-Lowerone,'',ll[1]-Lowertwo, '',ll[2]-Lowerone,'',ll[3]-Lowertwo)
           Dimrow=(ll[0]-Lowerone)*rangeP2+(ll[1]-Lowertwo)
           Dimcol=(ll[2]-Lowerone)*rangeP2+(ll[3]-Lowertwo)
           SMatrix=SparseMatrix(Msize, Msize, {(Dimrow,Dimcol):TermsCoeff[ii][1]})+SMatrix
    return SMatrix.rank()

 ## The function to calculation the SRV of the quantum state
def SchmidtRankVector(expr):  
    rho=toHH(expr)
    if rho==0 :return [0,0,0]
    else: return sorted([PartialTrace(rho,1),PartialTrace(rho,2),PartialTrace(rho,3)],reverse=True)
   
 

'''             
#good example: Reflection(OAMHolo(OAMHolo(OAMHolo(OAMHolo(LI(DCState,b,c),a,-2),a,-2),c,4),b,4),a)   trigger: FF1,[3,4]
#bad example:  OAMHolo(OAMHolo(BS(BS(OAMHolo(OAMHolo(OAMHolo(OAMHolo(DCState,c,-3),a,3),f,-5),c,2),a,c),a,b),b,-3),e,-4)  trigger: FF1,[1]
'''

# initialize experimental setup with SPDC of lower-OAM   
OAMorder=1  # OAM number from -OAMorder to OAMorder
crystal1=DownConvOAM(OAMorder,a,b)
crystal2=DownConvOAM(OAMorder,c,d)
DCState=expand((crystal1+crystal2)**2)
print('  wish you a nice day   ')
print('_________Small-SPDC: calculate the output state of a setup_________ ')
print('     ')

# this part takes very long for large setup and high OAM 
initial_state=expand(Reflection(OAMHolo(OAMHolo(OAMHolo(OAMHolo(LI(DCState,b,c),a,-2),a,-2),c,4),b,4),a)) 
#initial_state=expand(OAMHolo(OAMHolo(BS(BS(OAMHolo(OAMHolo(OAMHolo(OAMHolo(DCState,c,-3),a,3),f,-5),c,2),a,c),a,b),b,-3),e,-4))
print('*********** Post-selection  **************')

ReducedVV=Trigger(MakeFF2(initial_state),FF1,[3,4]) 
#ReducedVV=Trigger(MakeFF2(initial_state),FF1,[1]) # this also takes long for high OAM SPDC
print('Small_SPDC Output: ',ReducedVV)
print('     ')

ReducedVVList = TermsCoeffList(ReducedVV)
print('ReducedVVList ',ReducedVVList) 
print('     ') 

## store the FF2 FF3 FF4 in a list 
AllTypesOfFF2=AllTypesOfFFl(ReducedVVList,FF2,0)
AllTypesOfFF3=AllTypesOfFFl(ReducedVVList,FF3,1)
AllTypesOfFF4=AllTypesOfFFl(ReducedVVList,FF4,2)
print('AllTypesOfFF2:  ', AllTypesOfFF2,';')
print('AllTypesOfFF3:  ', AllTypesOfFF3,';')
print('AllTypesOfFF4:  ', AllTypesOfFF4,';')
print('     ')
DimVecEncode=sorted([len(AllTypesOfFF2),len(AllTypesOfFF3),len(AllTypesOfFF4)],reverse=True) ## get the maximum number of used dimensions

## check all the coefficients in the Output state are the same and the number of the term 
TermofState, CoeffsEqual=CheckTermCoeffOfState(ReducedVVList)
print('TermofState:  ',TermofState,' CoeffsEqual?  ', CoeffsEqual)
print('     ')

if DimVecEncode[0]>TermofState: ## check the number of term and dimensions
   InfoLetter1='_MoreTerms'
elif DimVecEncode[0]==TermofState:
   InfoLetter1='_OKTerms'   
    
if CoeffsEqual:  ## check whether all coefficients are the same
   InfoLetter2='_MaxEnt'
else:
   InfoLetter2='_NoMaxEnt'
   
InfoLetter=InfoLetter1+InfoLetter2
print('Small-SPDC Check Status: ',InfoLetter)
print('     ')

if DimVecEncode[2]>1 and InfoLetter=='_OKTerms_MaxEnt' :  #the index starts from 0
    DimVec=SchmidtRankVector(ReducedVV) 
    print('SRV: ',DimVec)
    print('     ')

if DimVec[2]>1 and TermofState== DimVec[0]:  ## check maximun entangled states
    OAMorder=5  # check high-order OAM modes
    crystal1=DownConvOAM(OAMorder,a,b)
    crystal2=DownConvOAM(OAMorder,c,d)
    DCState=expand((crystal1+crystal2)**2)
    
    print('_________Large-SPDC: calculate the output state of a setup_________ ')
    print('     ')
    # this part takes very long for large setup and high OAM 
    initial_state=expand(Reflection(OAMHolo(OAMHolo(OAMHolo(OAMHolo(LI(DCState,b,c),a,-2),a,-2),c,4),b,4),a))
#    initial_state=expand(OAMHolo(OAMHolo(BS(BS(OAMHolo(OAMHolo(OAMHolo(OAMHolo(DCState,c,-3),a,3),f,-5),c,2),a,c),a,b),b,-3),e,-4))

    print('*********** Post-selection  **************')
    
    ReducedVVHO=Trigger(MakeFF2(initial_state),FF1,[3,4])
#    ReducedVVHO=Trigger(MakeFF2(initial_state),FF1,[1])   # this also takes long for high OAM SPDC
    print('Large_SPDC Output (before checking):  ',ReducedVVHO)
    print('     ')
     
    ReducedVVHOList = TermsCoeffList(ReducedVVHO)
    print('ReducedVVHOList ',ReducedVVHOList) 
    print('     ') 
    
    ## store the FF2 FF3 FF4 in a list 
    AllTypesOfFF2Full=AllTypesOfFFl(ReducedVVHOList,FF2,0)
    AllTypesOfFF3Full=AllTypesOfFFl(ReducedVVHOList,FF3,1)
    AllTypesOfFF4Full=AllTypesOfFFl(ReducedVVHOList,FF4,2)
    print('AllTypesOfFF2Full:  ', AllTypesOfFF2Full,';')
    print('AllTypesOfFF2Full:  ', AllTypesOfFF3Full,';')
    print('AllTypesOfFF2Full:  ', AllTypesOfFF4Full,';')
    print('     ')
    
    ComplementDC2=list(AllTypesOfFF2Full.difference(AllTypesOfFF2))
    ComplementDC3=list(AllTypesOfFF3Full.difference(AllTypesOfFF3))
    ComplementDC4=list(AllTypesOfFF4Full.difference(AllTypesOfFF4))
    
    for ii in range(len(ComplementDC2)):
       ReducedVVHO=ReducedVVHO.replace(ComplementDC2[ii],0, map=False, simultaneous=True, exact=False)
    for ii in range(len(ComplementDC3)):
       ReducedVVHO=ReducedVVHO.replace(ComplementDC3[ii],0, map=False, simultaneous=True, exact=False)        
    for ii in range(len(ComplementDC4)):
       ReducedVVHO=ReducedVVHO.replace(ComplementDC4[ii],0, map=False, simultaneous=True, exact=False) 
       
    print('Large_SPDC Output (after checking): ',ReducedVVHO)
    print('   ')
    if ReducedVV==ReducedVVHO:
        print('This is a good SRV:', DimVec)
    else :
        print('This is a bad SRV beacuse of high-order OAM')
    
 
 
