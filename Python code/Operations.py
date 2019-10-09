# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 12:47:38 2019

@author: Miranda
"""
from sympy import mod_inverse as mminv
from sympy import isprime as isprime
import math
import random
import time
#Addition

def wadd(p1,p2,a,b,p):
    """Adds p1 and p2 on a WNF curve"""
    if math.isinf(p1[0]) and math.isinf(p2[0]):
        return math.inf,1 #O + O = O
    if math.isinf(p1[0]): return p2 #O + p2 = p2
    if math.isinf(p2[0]): return p1 #p1 + O = p1
    else:
        x1,y1 = p1 #assigns coordinates
        x2,y2 = p2 #assigns coordinates
        if y1**2 % p != (x1**3 + a*x1 + b) % p: #checks if point is on curve
            return False #otherwise reject
        if y2**2 % p != (x2**3 + a*x2 + b) % p: #checks if point is on curve
            return False #otherwise reject
        if x1==x2 and -y2 % p == y1:
            return math.inf,1 #(x1,y1) + (x1,-y1) = O
        elif x1==x2:
            #point doubling
            sl = ((3 * x1**2 + a)*mminv(2*y1,p)) % p #slope
            x = (sl**2 - x1 - x2) % p 
            return x, -(sl * (x - x1) + y1) % p
        else:
            #point addition
            sl = (y2-y1)*mminv(x2-x1,p) % p #slope
            x = (sl**2 - x1 - x2) % p
            return x, -(sl * (x - x1) + y1) % p
        
def madd(p1,p2,A,B,p):
    """Adds p1 and p2 on a Montgomery curve"""
    if B*(A**2 - 4) == 0 % p: #checks Montgomery requirements
        print("Singular Curve")
        return
    if math.isinf(p1[0]) and math.isinf(p2[0]):
        return math.inf,1 #O + O = O
    if math.isinf(p1[0]): return p2 #O + p2 = p2
    if math.isinf(p2[0]): return p1 #p1 + O = p1
    else:
        x1,y1 = p1 #assigns coordinates
        x2,y2 = p2 #assigns coordinates
        if B*y1**2 % p != (x1**3 + A*x1**2 + x1) % p: #checks if point is on curve
            return False #otherwise reject
        if B*y2**2 % p != (x2**3 + A*x2**2 + x2) % p: #checks if point is on curve
            return False #otherwise reject
        if x1==x2 and -y2 % p == y1:
            return math.inf,1 #(x1,y1) + (x1,-y1) = O
        elif x1==x2:
            #point doubling
            sl = ((3 * x1**2 + 2*A*x1 + 1)*mminv(2*B*y1,p)) % p #slope
            x = (B*(sl**2) - 2*x1 - A) % p
            return x, (sl*(x1 - x) - y1) % p
        else:
            #point addition
            sl = (y2-y1)*mminv(x2-x1,p) % p #slope
            x = (B*(sl**2) - x1 - x2 - A) % p
            return x, (sl*(x1 - x) - y1) % p
            
def eadd(p1,p2,a,d,p):
    """
    Adds p1 and p2 on a twsited Edwards curve
    We assume that there exists a pt of order 4, so identity is (0,1)
    """
    if d*(1 - d) == 0 % p: #checks Edwards requirements
        print("Singular Curve")
        return
    x1,y1 = p1
    x2,y2 = p2 #check if in curve here
    if (a*x1**2 + y1**2) % p != (1 + d*(x1**2)*(y1**2)) % p: #checks point is on the curve
        return False #otherwise reject
    if (a*x2**2 + y2**2) % p != (1 + d*(x2**2)*(y2**2)) % p: #checks point is on the curve
        return False #otherwise reject
    if p1 == (0,1) and p2 == (0,1):
        return 0,1 #O + O = O
    if p1 == (0,1): return p2 #O + p2 = p2
    if p2 == (0,1): return p1 #p1 + O = p2
    if y1==y2 and -x2 % p == x1: #(x1, y1) + (-x1, y1) = O
        return 0,1
    else:
        #point addition
        dxy = d*x1*x2*y1*y2 % p
        return (x1*y2 + y1*x2)*mminv(1+dxy,p) % p, (y1*y2 - a*x1*x2)*mminv(1-dxy,p) % p
    
#Multiplication
      
def wmul(p1,m,a,b,p):
    """Multiply p1 by m on a WNF curve"""
    """Works - prelim checks"""
    start = time.process_time() * 1000
    kp = math.inf,1 #starts at O
    q = p1 #starts at 1P
    k = bin(m)[2:] #converts m to binary
    for ii in reversed(k): #reverses binary digits
        if ii == '1': #if k(ii) = 1, add q to kp
            kp = wadd(kp,q,a,b,p)     
        q=wadd(q,q,a,b,p) #double q each iteration
    print(time.process_time()*1000 - start)
    return kp

def mmul(p1,m,A,B,p):
    """Multiply p1 by m on a Montgomery curve"""
    start = time.process_time() * 1000
    kp = math.inf,1 #starts at O
    q = p1 #starts at 1P
    k = bin(m)[2:] #converts m to binary
    for ii in reversed(k): #reverses binary digits
        if ii == '1': #if k(ii) = 1, add q to kp
            kp = madd(kp,q,A,B,p)     
        q=madd(q,q,A,B,p) #double q each iteration
    print(time.process_time() * 1000 - start)
    return kp

def emul(p1,m,a,d,p):
    """Multiply p1 by m on a twsited Edwards curve"""
    start = time.process_time() * 1000
    kp = 0,1 #starts at O
    q = p1 #starts at 1P
    k = bin(m)[2:] #converts m to binary
    for ii in reversed(k): #reverses binary digits
        if ii == '1': #if k(ii) = 1, add q to kp
            kp = eadd(kp,q,a,d,p)     
        q=eadd(q,q,a,d,p) #double q each iteration
    print(time.process_time()*1000 - start)
    return kp

#Key Production
    
def wkeys(g,n,a,b,p):
    d = random.randint(1,n-1) #n order of generator G
    H = wmul(g,d,a,b,p)
    return d, H

def mkeys(g,n,A,B,p):
    d = random.randint(1,n-1) #n order of generator G
    H = mmul(g,d,A,B,p)
    return d, H

def ekeys(g,n,a,d,p):
    d = random.randint(1,n-1) #n order of generator G
    H = emul(g,d,a,d,p)
    return d, H



