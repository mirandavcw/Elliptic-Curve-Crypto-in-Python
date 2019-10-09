# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:04:00 2019

@author: Miranda
"""

from sympy import mod_inverse as mminv
from sympy import isprime as isprime
import math
import random
import time
time.time() * 10000


def legendres(dis,p):
    """legendre symbol"""
    ls = pow(dis, (p-1)//2, p)
    if ls == -1 % p:
        return -1
    else: return ls 

print(legendres(12,13))

def modroot(a,p):
    """Finds the square root of a modulo p"""
    a = a%p #reduces mod p
    ls=legendres(a,p) #checks Legendre Symbol
    if isprime(p) == 0:
        print("p not prime") #rejects if p not prime
        return
    if ls != 1:
        print("No real root") #rejects if Legendre Symbol is -1
        return False
    if p % 4 == 3: #if p mod 4 = 3, then sqrt of a must be +/- a^(p+1)/4
        b = pow(a, (p+1)//4, p)
        return b % p, -b % p
    else: #otherwise, go to TS algorithm
        return tonellishanks(a,p)
    
def tonellishanks(a,p):
    """Tonelli-Shanks algorithm for modular square roots"""
    #function is only ever called by modroot, so we need not check the Legendre Symbol
    t = (p-1)//2 #initialises t
    s = 1 #initialises s
    while t % 2 == 0:
        t = t//2 #divides t by two until the odd factor is found
        s = s+1 #stores power of two
    u = 2 #initialises u
    while legendres(u,p) == 1:
        u = random.randint(0,p) #chooses new u until a non-square is found
    v = pow(u,t,p) #computes v
    b = pow(a,(t+1)//2, p) #computes first guess
    #Now we find m - one binary digit at a time
    m = 0 #initialises m
    for ii in range(s-1): #from 0 to s-2:
        mii = pow(a,t*2**(s-2-ii),p) #computes value to decide if m_ii is 1 or 0
        if mii == 1:    mi=0 #stores m_ii
        else:   mi=1
        m = (m + mi*(2**(ii))) % p #computes m up to most current binary digit
    b = b*pow(v,m,p) % p #computes positive solution
    if b**2 % p != a % p: #double check
        b=tonellishanks(a,p)[0] #repeat if incorrect
    return b, -b % p
    
def modcube(a,p):
    """Finds the cube root of a modulo p --- algorithm by Prof. Cameron"""
    if p % 3 == 1 and p % 9 != 1: #Cameron method
        c = pow(a, (p - 1)//3, p) #check if cube root exists
        if c != 1: return  False #if not, reject
        else:
            if p%9==4: #if p mod 9 = 4...
                k = (2*p + 1)//9 #k can be given as this
            elif p%9==7: #if p mod 9 = 7...
                k = (p+2)//9 #k can be given as this
            else:
                k = random.randint(0,(p+1)/2) #otherwise choose random k
                while True:
                    k = k+1 #increase k until...
                    if 3*k % ((p-1)//3) == 1: #suitable k is found
                        break
            b = pow(a,k,p) #compute cube root
        return b
    elif p % 3 == 2: #we can use short formula:
        b = pow(a, (2*p-1)//3, p)  
        return b
    
    
def w2check(a,b,p):
    """Checks for point of order two on WNF curve"""
    """Returns x coordinate alpha if exists"""
    dis = -16*(4 * a**3 + 27 * b**2) % p #calculates discriminant
    ldis = legendres(dis,p) #finds Legendre Symbol of discriminant
    if dis == 0:
        print("Discriminant = 0") #rejects if discriminant is 0
        return
    if ldis % 2 == 1:
        #1 root if -1, 3 if 1, doing both cases here
        r = ((b**2)*mminv(4,p) + (a**3)*mminv(27,p)) % p #using cubic formulae 
        rt = modroot(r,p) #sqrt of r
        if not rt:
            print("No transformation") #if sqrt does not exist, no alpha exists
            return
        else:
            rt = rt[0] #assign one root to use in algebra
        bf = -b*mminv(2,p) % p #cubic formulae
        alpha = (modcube(bf+rt,p) + modcube(bf-rt,p)) #cubic formulae
        if (alpha**3 + a*alpha + b) % p == 0:  
            return alpha % p #double checking value of alpha
        else:
            print("error") #reject otherwise, output alpha for error checking
        return alpha % p
    else:
        print("No point of order two")
        return 

#Points
    
def wtom(a,b,p):
    """Converts WNF curve to Montgomery"""
    """Checks for conversion"""
    start = time.process_time()*1000
    alpha = w2check(a,b,p) #finds point of order two on W
    if not alpha:
        print("No conversion") #if no such point, reject
        return
    else:
        if not modroot(3 * alpha**2 + a, p): #checks if square root exists
            print("No conversion") #otherwise, reject
            return
        B = mminv(max(modroot(3 * alpha**2 + a, p)), p) #calulates B
        A = 3*alpha*B % p #calulates A
        print(B,"y**2 = x**3 +",A,"x**2 + x") #prints curve equation
        print("An Isomorphism!") #all existing W to M mappings are isomorphisms
        print(time.process_time()*1000 - start)
        return A,B,alpha
   
def mtow(A,B,p):
    """Converts Montgomery curve to WNF"""
    #conversion always exists - no need to check
    start = time.process_time() *1000
    a = ((3-A**2)*mminv(3*B**2,p)) % p #calculates a
    b = ((2*A**3 - 9*A)*mminv(27*B**3,p)) % p #calculates b
    print("y**2 = x**3 + ",a, "x +", b) #prints curve equation
    print("An Isomorphism!") #all existing W to M mappings are isomorphisms
    print(time.process_time()*1000 - start)
    return a,b

def mtoe(A,B,p):
    """Converts Montgomery curve to twisted Edwards"""
    """Checks for conversion"""
    #start = time.process_time()*1000
    if A < 3 % p and A > -2 % p: #checks requirement of A 
        print("A in {-2,2}") #if fails, reject
        return
    elif B == 0 % p: #checks B non zero
        print("B = 0") #otherwise reject
        return
    a = ((A+2)*mminv(B,p)) % p #calculates a
    d = ((A-2)*mminv(B,p)) % p # calculates d
    print(a,"x**2 + y**2 = 1 + ",d,"x**2y**2") #prints curve equation
    if legendres(a,p) == legendres(d,p):
        print("Not an Isomorphism :(") #not isomorphic if a & d both square/not square
    elif legendres(a,p) == 1:
        print("An Isomorphism!") #isomorphic if a square, d not square
    elif legendres(a,p) == -1:
        print("Not an Isomorphism :(") #not isomorphic if a not square
   # print(time.process_time()*1000 - start)
    return a,d

def etom(a,d,p):
    """Converts twisted Edwards to Montgomery"""
    """Checks for conversion"""
    #start = time.process_time() *1000
    if a == d: #checks parameters distinct
        print("a=d") #otherwise reject
        return
    elif a == 0 or d == 0: #checks parameters non zero
        print("Must be non zero") #otherwise reject
        return
    A = (2*(a+d)*mminv(a-d,p)) % p #calculates A
    B = (4*mminv(a-d,p)) % p #calculates B
    print(B,"y**2 = x**3 +",A,"x**2 + x") #prints curve equation
    if legendres(a,p) == legendres(d,p):
        print("Not an Isomorphism :(") #not isomorphic if a & d both square/not square
    elif legendres(a,p) == 1: 
        print("An Isomorphism!") #isomorphic if a square, d not square
    #print(time.process_time()*1000 - start)
    return A,B

def etow(a,d,p):
    """Converts twisted Edwards to Weierstrass"""
    """Uses that W must have a pt of order 4 therefore T is transformable to M
    """
    A,B = etom(a,d,p)
    a,b = mtow(A,B,p)
    return a,b

def wtoe(a,b,p):
    """Converts Weierstrass to twisted Edwards """
    """Uses that W must have a pt of order 4 therefore T is transformable to M
    """
    start = time.process_time() *1000
    A,B,al = wtom(a,b,p)
    a,d = mtoe(A,B,p)
    print(time.process_time()*1000 - start)
    return a,d

#Curves
def pwtom(p1,a,b,p,A,B):
    """Transforms a Weierstrass point to a Montgomery point"""
    x1,y1 = p1 #assigns coordinates
    if math.isinf(p1[0]): return math.inf,1 #inf -> inf
    if y1**2 % p != (x1**3 + a*x1 + b) % p: #checks if point is on curve
        return False #otherwise reject
    if not A: #if user does not input A, calulates it
        s,alpha = wtom(a,b,p)[1],wtom(a,b,p)[2] #finds B, alpha from wtom
        if not alpha: #checks if a point of order two exists
            print("No conversion") #otherwise reject
            return
    else:
        alpha = A*mminv(3*B,p) % p #if user inputs A,B calculate alpha from them
        s = B #as above
    #exceptional points
    if p1 == (alpha,0): return 0,0 #if point has order two, map to 0,0
    #general points
    else:
        return s*(x1-alpha) % p, s*y1 % p 
    
def pmtow(p1,A,B,p):
    """Transforms a Montgomery point to a Weierstrass point"""
    start = time.process_time()*1000
    x1,y1 = p1 #assigns coordinates
    if B*y1**2 % p != (x1**3 + A*x1**2 + x1) % p: #checks if point is on curve
        return False #otherwise reject
    #exceptional points
    if math.isinf(x1): 
        print(time.process_time()*1000 - start)
        return math.inf,1 #inf -> inf
    #general points
    else:
        print(time.process_time()*1000 - start)
        return (x1*mminv(B,p) + A*mminv(3*B,p)) % p, y1*mminv(B,p) % p
    
def pmtoe(p1,A,B,p):
    """Transforms a Montgomery point to a twisted Edwards point"""
    start = time.process_time()*1000
    if A < 3 % p and A > -2 % p: #checks requirements on A for curve mapping
        print("A in {-2,2}") #if fails, reject
        return False
    elif B == 0 % p: #checks requirements on B for curve mapping
        print("B = 0") #if fails, reject
        return False
    x1,y1 = p1 #assigns coordinates
    if math.isinf(x1): return 0,1 #inf -> (0,1)    
    if B*y1**2 % p != (x1**3 + A*x1**2 + x1) % p: #checks if point is on curve
        return False #otherwise reject
    #exceptional points
    elif y1 == 0 % p: 
        print(y1)
        if x1 == 0 % p: return 0, -1 % p #(0,0) -> (0,-1)
        elif 2*x1 % p == (-A + modroot((A+2)*(A-2),p)[0]) % p:
                return math.inf,2 #point of order two on M -> inf ('positive' root)
        elif 2*x1 % p == (-A - modroot((A+2)*(A-2),p)[0]) % p:
                return math.inf,-2 #point of order two on M -> inf ('negative' root)
        else: return
    elif x1 == -1 % p: 
        if y1 == modroot((A-2)*mminv(B,p),p):
            return math.inf,4 #point of order four on M -> inf ('positive' root)
        elif y1 == -modroot((A-2)*mminv(B,p),p) % p:
            return math.inf,-4 #point of order four on M -> inf ('negative' root)
    #general points
    else:
        print(time.process_time()*1000 - start)
        return x1*mminv(y1,p) % p, (x1-1)*mminv(x1+1,p) % p

def petom(p1,a,d,p):
    """Transforms a twisted Edwards point to a Montgomery point"""
    if a == d: #checks curve mapping requirements
        print("a=d") #if fails, reject
        return False
    elif a == 0 or d == 0:
        print("Must be non zero")
        return False
    x1,y1 = p1 #assigns coordinates
    #exceptional points
    if math.isinf(x1):
        A = (2*(a+d)*mminv(a-d,p)) % p #calculates A
        B = (4*mminv(a-d,p)) % p #calculates B
        if y1 % p == 2:
            x = (-A + modroot((A+2)*(A-2),p)[0])*mminv(2,p) % p #point of order two on M ('positive' root)
            return x,0
        elif y1 % p == -2 % p:
            x = (-A - modroot((A+2)*(A-2),p)[0])*mminv(2,p) % p #point of order two on M ('negative' root)
            return x,0
        elif y1 % p == 4:
            y = modroot((A-2)*mminv(B,p),p) #point of order four on M ('positive' root)
            return -1,y
        elif y1 % p == -4 % p:
            y = -modroot((A-2)*mminv(B,p),p) % p #point of order four on M ('negative' root)
            return -1,y
        else: return #no other points at infinity should exist on a twisted Edwards curve
    if (a*x1**2 + y1**2) % p != (1 + d*(x1**2)*(y1**2)) % p: #checks point is on the curve
        return False #otherwise reject
    elif x1 % p == 0 and y1 % p == 1: return math.inf,1 #(0,1) -> inf
    elif x1 % p == 0 and y1 % p == -1 % p: return 0,0 #(0,-1) -> (0,0)
    #general points
    else:
        return (1+y1)*mminv(1-y1,p) % p, (1+y1)*mminv(x1*(1-y1),p) % p
        
def petow(p1,a,d,p):
    """Transforms a twisted Edwards point to a Weierstrass point"""
    if a == d: #checks curve mapping requirements
        print("a=d") #if fails, reject
        return
    elif a == 0 or d == 0:
        print("Must be non zero")
        return
    x1,y1 = p1 #assigns coordinates
    #exceptional points
    if math.isinf(x1):
        #computes inf seperately, since may be (inf,4) or (inf,2)
        p2 = petom(p1,a,d,p) #computes point at infinity -> point on montgomery
        A,B = (2*(a+d)*mminv(a-d,p)) % p, (4*mminv(a-d,p)) % p #calculates A and B
        p3 = pmtow(p2,A,B,p) #computes point on montgomery -> point on W
        return p3
    if (a*x1**2 + y1**2) % p != (1 + d*(x1**2)*(y1**2)) % p: #checks point is on curve
        return False #otherwise reject
    elif x1 % p == 0 and y1 % p == 1:
        return math.inf,1 #(0,1) -> inf -> inf
    else:
        p2 = petom(p1,a,d,p) #computes point on Montgomery
        A,B = (2*(a+d)*mminv(a-d,p)) % p, (4*mminv(a-d,p)) % p #calculates A and B
        p3 = pmtow(p2,A,B,p) #computes point on W
        return p3
    
def pwtoe(p1,a,b,p):
    """Transforms a Weierstrass point to a Montgomery point"""
    """Need to check req for W (pt of order 4 exists etc)
    """
    x1,y1 = p1 #assigns cordinates
    if math.isinf(x1):
        return 0,1 #inf -> inf -> (0,1)    
    if y1**2 % p != (x1**3 + a*x1 + b) % p: #checks point is on curve
        return False #otherwise reject
    else:
        A,B,al = wtom(a,b,p) #calculate Montgomery parameters
        p2 = pwtom(p1,a,b,p,A,B) #calculate point on Montgomery
        p3 = pmtoe(p2,A,B,p) #calculate point on twisted Edwards
        return p3

        
       
            
            
    
