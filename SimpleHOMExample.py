# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:23:21 2018

@author: MG
"""

 
import sympy as sp
 
# path information 
#a = sp.IndexedBase("a")
#b = sp.IndexedBase("b")
#c = sp.IndexedBase("c")
#d = sp.IndexedBase("d")

a, b, c, d, e, f= sp.symbols('a b c d e f', cls=sp.IndexedBase)
FF1, FF2, FF3, FF4, FF5, FF6= sp.symbols('FF1 FF2 FF3 FF4 FF5 FF6', cls=sp.IndexedBase)
x1, x2, x3, x4, coff= sp.symbols('x1 x2 x3 x4 coff', cls=sp.Wild)

#polarization
H,V=sp.symbols('H V')
Modelist=[H,V] 


"""
Polarization
simpleHOMExample

-> Hong-Ou-Mandel destructive interference
-> Hong-Ou-Mandel constructive interference
-> 4-photon GHZ generation in postselection
 
"""
def DownConvPol(p1,p2):
    expr=p1[H]*p2[H] + p1[V]*p2[V]
    return expr

def HWP(expr,p):
    #expr=expr.subs([(p[H],p[V]),(p[V],p[H])],simultaneous=True) 
    expr=expr.subs([(p[H],p[V]),(p[V],-p[H])],simultaneous=True) # in your code, V change to -H
    return expr

def BSPol(expr,p1,p2): 
    expr=expr.subs([(p1[x1],(p2[x1]+sp.I*p1[x1])/sp.sqrt(2)) for x1 in Modelist]+[(p2[x1],(p1[x1]+sp.I*p2[x1])/sp.sqrt(2)) for x1 in Modelist],simultaneous=True)
    return expr
#
def PBS(expr,p1,p2):
    expr=expr.subs([(p1[H],p2[H]),(p1[V],sp.I*p1[V]),(p2[H],p1[H]),(p2[V],sp.I*p2[V])],simultaneous=True)
    return expr

# replaceRule --- _ Blank() and replace rules in mathematica; 
# cannot be used to BS and othert because of the simultaneous=True, need to be replace in once.
def replaceRule(expr, repls):
    for k, m in repls.items():
        expr = expr.replace(k, m, map=False, simultaneous=True, exact=False)   
    return expr

def MakeFF(expr): 
    NFoldrepls = {coff*a[x1]*b[x2]*c[x3]*d[x4] : coff*FF1[x1]*FF2[x2]*FF3[x3]*FF4[x4], coff*a[x1]*a[x2] : 0, coff*b[x1]*b[x2] : 0, coff*c[x1]*c[x2] : 0, coff*d[x1]*d[x2] : 0}
    expr=replaceRule(expr,NFoldrepls)
    return expr

 
"""
We use a SPDC crystal which creates an entangled |phi+> state. 
Afterwards, we put a beam splitter. We will see Hong-Ou-Mandel interference 
(i.e. all photons in one the same output arm):
"""
print('Hong-Ou-Mandel destructive interference') 
print('                      ') 
psi=DownConvPol(b,c)
print('initial state:  psi=',psi)
 

psi2=sp.expand(BSPol(psi,b,c))
print('after BS:  psi2=',psi2)
 

"""
Now we again create a SPDC, and change put a polariser at 90\[Degree] in arm b, 
and add a phase of i in arm b, leading to a |psi-> state:
"""  
print('                      ') 
print('Hong-Ou-Mandel constructive interference') 
print('                      ') 
print('**********function calls step by step************') 
psi=DownConvPol(a,b)
print('initial state:  psi=',psi) 

psi2=HWP(psi,a)
print('after the HWP in path a:  psi2=',psi2)

psi3=sp.expand(BSPol(psi2,a,b))
print('psi2 after the BS:  psi3=',psi3)
 

#the same can be achieved in a sequence of function calls
psi3 = sp.expand(BSPol(HWP(DownConvPol(a,b),a),a,b))
print('                      ') 
print('********a sequence of function calls**********')
print('with HWP in path a and then BS, psi3=',psi3)
 


"""
The final example shows how the 4-photon Polarisation GHZ generation 
(using a PBS) can be described:
"""
print('                      ') 
print('4-photon GHZ generation in postselection')
print('                      ') 
 
psi=DownConvPol(a,b)*DownConvPol(c,d)
print('initial state:  psi=',psi)
 
psi2=sp.expand(PBS(psi,b,c)) 
print('after BS for path b and c:  psi2=',psi2)

psi3=MakeFF(psi2)
print('4-fold coincidence:  psi3=',psi3)
 



