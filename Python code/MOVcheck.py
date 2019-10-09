# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 14:07:24 2019

@author: Miranda
"""

def MOVcheck(n,p):
    for k in range(1,20):
        q = (p**k - 1) % n
        if q == 0:
            print(k)
            return "oh no!"
    return "Safe for all k < 20"
        
def findaG(n,a,b,p):
    while True:
        possx = random.randint(0,n-1)
        possy2 = (pow(possx,3,p) + possx*a + b) % p
        possy=modroot(possy2,p)
        if not possy:
            possy=1
        else:
            possG=possx,possy[0]
            G = wmul(possG,n,a,b,p)
            if math.isinf(G[0]):
                break
    return possG
            
def whatisorder(a,b,p):
    while True:
        possx = random.randint(0,p-1)
        possy2 = (pow(possx,3,p) + possx*a + b) % p
        possy=modroot(possy2,p)
        if not possy:
            possy=1
        else:
            possG=possx,possy[0]
            break
    P = possG
    n=1
    while True:
        P=wadd(possG,P,a,b,p)
        n=n+1
        if math.isinf(P[0]):
            break
    return possG, n
            
def timethelads(n,p):
    tot=0
    for ii in range(10):
        end=rootytooty(n,p)
        tot=tot+end
    return tot/10

def howmanyones(a):
    b=bin(a)
    c=0
    for ii in b:
        if ii == '1':
            c=c+1
    return c

def rootytooty(n,p):
    x=()
    while len(x)<n:
        r=random.randint(1,p-1)
        if legendres(r,p) == 1:
            x=x+(r,)
    start = time.process_time() * 1000
    for ii in range(len(x)):
        modroot(x[ii],p)
    end = time.process_time()*1000 - start
    return end