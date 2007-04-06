import sys

import numarray

import chemev
from chemev import iso,optimization,utils

import time

def fitting_iteration(pars,iter=10):
    bestfit = [None]
    ndata=sum(data.flat)
    def f(params):
        m=utils.normalize(iso.computeCMD(params,isos),ndata)
        return utils.loglikelihood(m,data)

    d = numarray.maximum(data,1e-20)
    llhC=sum( (d*numarray.log(d)).flat )
    def b(params,value,iter):
        "Print status everytime a best fit is found."
        m=iso.computeCMD(params,isos)
        C=sum(data.flat)/sum(m.flat)
        x=numarray.array(params)
        bestfit[0] = x*C
        file("bestfit117","w").write(str(list(bestfit[0])))
        print "bestfit in the file is:",f(bestfit[0])
        print "henry:",value,"tom:",2.0*(value+llhC),"iter:",iter,"norm:",C

    #optimization.minc(optimization.fmin_scipy_l_bfgs_b,f,pars,callback=b,iter=iter)
    x,fv,g= optimization.fmin_scipy_l_bfgs_b(f,pars,iter=iter)
    b(x,fv,-1)
    return x

print "reading isochrones"
data=iso.readfits(chemev.isodir+"/117/halo/datarr.fits")
isos = iso.readisos(chemev.isodir+"/117/halo")

try:
    pars=eval(open("bestfit117").readline())
except IOError:
    pars=[1.0]*117

print "fitting"

#print pars
pars=fitting_iteration(pars,iter=50000/119)
#pars=fitting_iteration(pars,iter=10)
#print pars

print "done"
