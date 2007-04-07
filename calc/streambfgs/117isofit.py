import sys

#import numarray

import chemev
from chemev import iso,optimization,utils

from numpy import array, maximum, log

import time

def fitting_iteration(pars,iter=10):
    bestfit = [None]
    ndata=sum(data.flat)
    def f(params):
        m=utils.normalize(iso.computeCMD(params,isos),ndata)
        return utils.loglikelihood(m,data)

    d = maximum(data,1e-20)
    llhC=sum( (d*log(d)).flat )
    def b(params,value,iter):
        "Print status everytime a best fit is found."
        m=iso.computeCMD(params,isos)
        C=sum(data.flat)/sum(m.flat)
        x=array(params)
        bestfit[0] = x*C
        file("bestfit117","w").write(str(list(bestfit[0])))
        print "henry:",value,"tom:",2.0*(value+llhC),"iter:",iter,"norm:",C

    a= optimization.minc(optimization.fmin_bfgs,f,pars,callback=b,iter=iter)
    b(a, f(a), -1)
    return a


print "reading isochrones"
data=iso.readfits(chemev.isodir+"/117/stream/datarr.fits")
isos = iso.readisos(chemev.isodir+"/117/stream")

pars=[1.0]*117

print "fitting"

pars=fitting_iteration(pars,200)

print "done"
