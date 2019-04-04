# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 23:36:08 2019

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

a, b, c, d, e, f, FF1, FF2, FF3, FF4, FFn, HH, GG1, GG2, GG3, GG4=map(sp.IndexedBase,['a','b','c','d','e', 'f', 'FF1','FF2','FF3','FF4','FFn','HH','GG1','GG2','GG3','GG4'])  
l,l1, l2, l3, l4, l5, l6, l7, l8, x1, x2, x3, x4, x5, x6, coeff, powern =map(sp.Wild,['l','l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'coeff', 'powern'])

zero=Symbol('zero') ## using in 4-fold coincidence 
psi=Symbol('psi')  ## represent a quantum state

## corresponding symbolic number 
sqr2=sqrt(2)/2  ## equals to 1/sqrt(2)
imagI=I ## equals to Imaginary
Pi=pi

PossiblePath=[('a','b'),('a','c'),('a','d'),('a','e'),('a','f'),('b','c'),('b','d'),('b','e'),('b','f'),('c','d'),('c','e'),('c','f'),('d','e'),('d','f'),('e','f')]
PossiblePathNum=['a','b','c','d','e','f'] 
 
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
    expr=expr.replace(p[l1],imagI*p[-l1], map=False, simultaneous=True, exact=False)
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

def Allsubsets(iterable): ## create whole strings from device number
    xs = list(iterable)
    return list(chain.from_iterable(combinations(xs,n) for n in range(len(xs)+1)))


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

def AllTypesOfFFn(expr,FFn):   ## find all the used mode in path (no repeated modes)
    dictadd=collect(expr, [FFn[x1]], evaluate=False)
    TermsCoeff=list(dictadd.items())
    NumOfPath=[]
    for ii in range(len(TermsCoeff)):
        HHList=TermsCoeff[ii][0]
        NumOfPath.append(HHList.indices[0])
    return set(NumOfPath)  ## returen used modes in patha

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
   
  
## for creating random setup
def DefineActions():  ## creat the possible actions
    actions=[]
    for ii in range(len(PossiblePath)):
        PosA = PossiblePath[ii][0]
        PosB = PossiblePath[ii][1] 
        actions.append("BS(XXX,"+PosA+","+PosB+")")
        actions.append("LI(XXX,"+PosA+","+PosB+")")
        
    for ii in range(6):
        Pos=PossiblePathNum[ii]
        actions.append("Reflection(XXX,"+Pos+")")
        actions.append("DP(XXX,"+Pos+")")    
        nHOM=5
        HOM_list=list(range(-nHOM,nHOM+1))
        HOM_list.remove(0)  ## [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
        for ii in range(1,len(HOM_list)+1):  ## index strart from 1
             actions.append("OAMHolo(XXX,"+Pos+","+str(ii)+")") 
    return actions


def SetupStrList(actionTable,actionlist):## returen setup list
    rvATL=[]
    for iii in range(len(actionlist)):
      rvATL.append(actionTable[actionlist[iii]]) 
    YYY='XXX'
    for kk in range(len(rvATL)-1,-1,-1): # in the inverted sequence
        YYY=YYY.replace('XXX',rvATL[kk])
    return YYY

def PrintOptions(expr,printflag): ## print   
    if printflag==1:
        print(expr)


PrintFlag=1 ## chose to give print out and only give the correct SRV

actions=DefineActions()
random.seed()
SizeOfAllList=1000000
LoopNum=1000000000
for iiLoop in range(LoopNum):  ## create random setup
    iiNum=iiLoop%10 + 6
    AllListNum=choice(range(len(actions)), size=(SizeOfAllList,iiNum))
    PrintOptions('**************************************************',PrintFlag)
    PrintOptions('-------create random lists of experiments-------',PrintFlag)
    PrintOptions('     ',PrintFlag)  
    
    for iiLists in range(len(AllListNum)):# try every setuplist
       PrintOptions('-------start to try different experiments-------',PrintFlag)
       PrintOptions('     ',PrintFlag)
        
       setupStr=SetupStrList(actions, AllListNum[iiLists])  ## get the setup string
       
       OAMorder=1  # OAM number from -OAMorder to OAMorder
       crystal1=DownConvOAM(OAMorder,a,b)
       crystal2=DownConvOAM(OAMorder,c,d)
       DCState=(crystal1+crystal2)**2  ## produce photon pairs form SPDC

       FuncStr=setupStr.replace("XXX",str(DCState)) ## give the output state of the created setup
       output_state=expand(eval(FuncStr)) ## you can also use sp.sympify(FuncStr,mylist), don't know which is fast
       outputState=MakeFF2(output_state) 
      
       AllTypesOfFF1=AllTypesOfFFn(outputState,FF1)
       AllCombfFF1=Allsubsets(AllTypesOfFF1)  ## get all the subsets of FF1 for trigger

       
       for ii in range(1,len(AllCombfFF1)):  ## try every trigger
           PrintOptions('-------start to try different triggers in Path A (FF1)-------',PrintFlag)
           PrintOptions('     ',PrintFlag)
           
           lowoutput_state=outputState
           ReducedVV=Trigger(lowoutput_state,FF1,list(AllCombfFF1[ii])) 
           PrintOptions('Small_SPDC Output State for the setup: '+str(ReducedVV),PrintFlag)
           PrintOptions('     ',PrintFlag)     
           
           AllTypesOfFF2=AllTypesOfFFn(ReducedVV,FF2) ## store the FF2 FF3 FF4 in a list 
           AllTypesOfFF3=AllTypesOfFFn(ReducedVV,FF3)
           AllTypesOfFF4=AllTypesOfFFn(ReducedVV,FF4)
        
           PrintOptions('AllTypesOfFF2: FF2 '+str(AllTypesOfFF2)+';',PrintFlag)
           PrintOptions('AllTypesOfFF3: FF3 '+str(AllTypesOfFF3)+';',PrintFlag)
           PrintOptions('AllTypesOfFF4: FF4 '+str(AllTypesOfFF4)+';',PrintFlag)
           PrintOptions('     ',PrintFlag)
           DimVecEncode=sorted([len(AllTypesOfFF2),len(AllTypesOfFF3),len(AllTypesOfFF4)],reverse=True) ## get the maximum number of used dimensions
        
          ## check all the coefficients in the Output state are the same and the number of the term 
           ReducedVVList = TermsCoeffList(ReducedVV)
           TermofState, CoeffsEqual=CheckTermCoeffOfState(ReducedVVList)
           PrintOptions('TermofState:  '+str(TermofState)+' CoeffsEqual (0: No; 1: Yes)? '+str(CoeffsEqual),PrintFlag)
           PrintOptions('     ',PrintFlag)
           
           if DimVecEncode[0]<TermofState: ## check the number of term and dimensions 
               InfoLetter1='_MoreTerms'
           elif DimVecEncode[0]==TermofState:
               InfoLetter1='_OKTerms'  
            
           if CoeffsEqual:  ## check whether all coefficients are the same
               InfoLetter2='_MaxEnt'
           else:
               InfoLetter2='_NoMaxEnt'
           
           InfoLetter=InfoLetter1+InfoLetter2
           PrintOptions('Small-SPDC Check Status: '+InfoLetter+'; DimVecEncode:'+str(DimVecEncode),PrintFlag)
           PrintOptions('     ',PrintFlag)
        
           if DimVecEncode[2]>1 and InfoLetter=='_OKTerms_MaxEnt' :  #the index starts from 0
               DimVec=SchmidtRankVector(ReducedVV) 
               PrintOptions('SRV: '+str(DimVec),PrintFlag)
               PrintOptions('     ',PrintFlag)
               PrintOptions('The corresponing setup: '+setupStr+'; FF1: '+str(list(AllCombfFF1[ii])),PrintFlag)
               PrintOptions('     ',PrintFlag)
               
               if DimVec[2]>1 and TermofState== DimVec[0]:  ## check maximun entangled states
                   
                   PrintOptions('No problems in lower OAM orders, Now let us check higher order OAM SPDC...',PrintFlag)
                   PrintOptions('     ',PrintFlag)
                   OAMorder=5  # check high-order OAM modes
                   crystal1=DownConvOAM(OAMorder,a,b)
                   crystal2=DownConvOAM(OAMorder,c,d)
                   DCState=(crystal1+crystal2)**2
                
                   FuncStr=setupStr.replace("XXX",str(DCState)) ## give the output state of the created setup
                   highoutput_state=expand(eval(FuncStr))
                   highoutputState=MakeFF2(highoutput_state) 
                   
                   ReducedVVHO=Trigger(highoutputState,FF1,list(AllCombfFF1[ii]))
                 
                   AllTypesOfFF2Full=AllTypesOfFFn(ReducedVVHO,FF2) ## store the FF2 FF3 FF4 in a list 
                   AllTypesOfFF3Full=AllTypesOfFFn(ReducedVVHO,FF3)
                   AllTypesOfFF4Full=AllTypesOfFFn(ReducedVVHO,FF4)        
                   PrintOptions('AllTypesOfFF2Full: FF2 '+str(AllTypesOfFF2Full)+';',PrintFlag)
                   PrintOptions('AllTypesOfFF3Full: FF3 '+str(AllTypesOfFF3Full)+';',PrintFlag)
                   PrintOptions('AllTypesOfFF4Full: FF4 '+str(AllTypesOfFF4Full)+';',PrintFlag)
                   PrintOptions('     ',PrintFlag)
                
                   ComplementDC2=list(AllTypesOfFF2Full.difference(AllTypesOfFF2))
                   ComplementDC3=list(AllTypesOfFF3Full.difference(AllTypesOfFF3))
                   ComplementDC4=list(AllTypesOfFF4Full.difference(AllTypesOfFF4))
                   if ComplementDC2 != []:
                       for iii in range(len(ComplementDC2)):
                           ReducedVVHO=ReducedVVHO.replace(FF2[ComplementDC2[iii]],0, map=False, simultaneous=True, exact=False)
                   if ComplementDC3 != []:
                       for iii in range(len(ComplementDC3)):
                           ReducedVVHO=ReducedVVHO.replace(FF3[ComplementDC3[iii]],0, map=False, simultaneous=True, exact=False)  
                   if ComplementDC4 != []:
                       for iii in range(len(ComplementDC4)):
                           ReducedVVHO=ReducedVVHO.replace(FF4[ComplementDC4[iii]],0, map=False, simultaneous=True, exact=False) 
                   
                   PrintOptions('Large_SPDC Output State for the setup: '+str(ReducedVVHO),PrintFlag)
                   PrintOptions('   ',PrintFlag)
                   if ReducedVV==ReducedVVHO:
                       print('This is a good SRV:', DimVec)
                       inforStr= 'SRV:'+str(DimVec)+'; Setup: '+setupStr+'; FF1:'+str(list(AllCombfFF1[ii]))+'\r\n'
                       print('info: ',inforStr)
                       goodlist = open("Goodcases.txt", "a+")
                       goodlist.write(inforStr)
                       goodlist.close()
                   else :
                       print('This is a bad SRV beacuse of high-order OAM', DimVec)
                       inforStr= 'SRV:'+str(DimVec)+'; Setup: '+setupStr+'; FF1:'+str(list(AllCombfFF1[ii]))+'\r\n'
                       print('info: ',inforStr)
                       badlist = open("Badcases.txt", "a+")
                       badlist.write(inforStr)
                       badlist.close()                                 
