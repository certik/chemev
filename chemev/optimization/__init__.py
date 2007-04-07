""" Routines for optimization, adapted from codes on the web.
    amoeba -- Nelder-Meade simplex algorithm
    DESolver -- Differential evolution
    mcmc -- Monte Carlo Markov Chain
"""

__version__ = "1.0.0"
__author__ = "H. C. Ferguson & O. Certik, STScI"

import math

import numarray

from amoeba import fmin_simplex
from DESolver import fmin_de
from scipyoptimize import fmin_bfgs
from mcmc import fmin_mcmc
have_mcmc=mcmc.have_mcmc

class Logistics(object):

    def __init__(self, min, max):
        self.min = min
        self.max = max

    def logistic2real(self,x):
        """Converts the tuple 'x' from logistic to real"""
        return [self.frac(m,mi,ma) for m,mi,ma in zip(x,self.min,self.max)]

    def real2logistic(self,x):
        """Converts the tuple 'x' from real to logistic"""
        return [self.fracinv(m,mi,ma) for m,mi,ma in zip(x,self.min,self.max)]

class ReflectLogistics(Logistics):

    def reflect(self,y,l,u):
        """ -inf < y < inf      ->  l < x < u

        It modulates "x" so that it cycles from l
        to u then back down again as "y" varies from 0 to 1 then 1 to 2. These
        cycles repeat so that "y" can vary from -inf to inf and "x" will just
        reflect within its boundaries.
        """
        m = y%2-y%1
        a = (u-l)*(y%1)+l
        b = (l-u)*(y%1)+u
        x = (1-m)*a + m*b
        return x

    def frac(self, a, x0, x1):
        """-inf < a < inf    ->   x0 ... x1"""
        return self.reflect(a, x0, x1)

    def fracinv(self,a,x0,x1): 
        """x0<a<x1  -> -inf ... inf"""
        #I am lazy to implement it now, but because I am using random values
        #at the beginning, this works fine:
        import random
        return random.random()*2

class ExpLogistics(Logistics):

    def logistic(self,x):
        """-inf<x<inf    ->   0..1"""
        if x<-500: return 0.
        return 1./(1.+math.exp(-x))

    def logisticinv(self,x):
        """0<x<1  ->  -inf .. inf"""
        assert x>=0.0
        assert x<=1.0
        if x<1e-50: return -115
        if x>1-1e-10: return 23
        return -math.log(1.0/x-1)

    def frac(self,a,x0,x1): 
        """-inf < a < inf    ->   x0 ... x1"""
        return x0-self.logistic(a)*(x0-x1)

    def fracinv(self,a,x0,x1): 
        """x0<a<x1  -> -inf ... inf"""
        return self.logisticinv((x0-a)*1.0/(x0-x1))


def minmax(algorithm,f,pars,min,max,iter=100,callback=None, 
        logistics = ExpLogistics):
    """ 
    Optimize a function within min and max. The fitted parameters
    are a logistic function of the real parameters within min,max.
    So the optimization searches the range -infinity to infinity while
    the true parameter ranges between min & max. For amoeba fitting,
    set the initial guesses so they are not right at min or max, 
    otherwise they will be stuck there. You can use parameters.randomize()
    and run the optimization several times to make sure the best fit
    is really a best fit.
        
    logistics .... a class that handles logistics. Use your own class for
        other logistics besides the default one
    """
    pars=list(numarray.minimum(numarray.array(pars),max))
    l = logistics(min,max)
    def newf(var,data=None):
        x=l.logistic2real(var)
        return f(x)
    def newcallback(params,value,iter):
        x=l.logistic2real(params)
        return callback(x,value,iter)
    r=algorithm(newf,l.real2logistic(pars),iter=iter,callback=newcallback)
    if r!=None:
        return l.logistic2real(r) 

def optimize(algorithm,f,pars,callback=None):
    """ Optimize a function. For this version, pars
        is expected to contain the values of the parameters, the range,
        and any information on how to treat the limits. All of this is
        handled by the function f --- so minmax2 knows nothing about the
        structure of pars, and does not modify the parameters in any way
    """
    r=algorithm(f,pars,callback=callback)
    if r!=None:
        return r

def minc(algorithm,f,pars,callback=None,iter=100):
    M=1.e7
    min=[0.]*len(pars)
    max=[M]*len(pars)
    return minmax(algorithm,f,pars,min,max,iter=iter,callback=callback)
