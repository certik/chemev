import sys

import numarray

import chemev
from chemev import iso,optimization,utils

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
        print "henry:",value,"tom:",2.0*(value+llhC),"iter:",iter,"norm:",C

    optimization.minc(optimization.fmin_simplex,f,pars,callback=b,iter=iter)
    return bestfit[0]

print "reading isochrones"
data=iso.readfits(chemev.isodir+"/117/halo/datarr.fits")
isos = iso.readisos(chemev.isodir+"/117/halo")

pars=[1.0]*117

print "fitting"

#this is the hardwired empirical key to achieve the best fit. :)
best_fit_path = [1000, 20, 30]
for iter in best_fit_path:
    print "### doing %d iterations ###"%iter
    pars=fitting_iteration(pars,iter)

print "done"
