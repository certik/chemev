import math

import numarray

from amoeba import fmin_simplex
from DESolver import fmin_de

def logistic(x):
    """-inf<x<inf    ->   0..1"""
    if x<-500: return 0.
    return 1./(1.+math.exp(-x))

def logisticinv(x):
    """0<x<1  ->  -inf .. inf"""
    assert x>=0.0
    assert x<=1.0
    if x<1e-50: return -115
    if x>1-1e-10: return 23
    return -math.log(1.0/x-1)

def frac(a,x0,x1): 
    """-inf < a < inf    ->   x0 ... x1"""
    return x0-logistic(a)*(x0-x1)

def fracinv(a,x0,x1): 
    """x0<a<x1  -> -inf ... inf"""
    return logisticinv((x0-a)*1.0/(x0-x1))

def logistic2real(x,min,max):
    """Converts the tuple 'x' from logistic to real"""
    return [frac(m,mi,ma) for m,mi,ma in zip(x,min,max)]

def real2logistic(x,min,max):
    """Converts the tuple 'x' from real to logistic"""
    return [fracinv(m,mi,ma) for m,mi,ma in zip(x,min,max)]

def minmax(algorithm,f,pars,min,max,callback=None):
    pars=list(numarray.minimum(numarray.array(pars),max))
    def newf(var,data=None):
        x=logistic2real(var,min,max)
        return f(x)
    def newcallback(params,value,iter):
        x=logistic2real(params,min,max)
        return callback(x,value,iter)
    r=algorithm(newf,real2logistic(pars,min,max),callback=newcallback)
    return logistic2real(r,min,max) 

def minc(algorithm,f,pars,callback=None):
    M=1.e7
    min=[0.]*len(pars)
    max=[M]*len(pars)
    return minmax(algorithm,f,pars,min,max,callback)
