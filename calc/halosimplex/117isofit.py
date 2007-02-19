import sys

import numarray

from chemev import iso,optimization,utils

print "reading isochrones"
#isodir="m31iso"
#isodir="m117halo"
#isodir="/home/ondrej/data/isochrones/117/halo"
isodir="../../isochrones/117/halo"
data=iso.readfits(isodir+"/datarr.fits")
isos = iso.readisos(isodir)

try:
    pars=eval(open("bestfit117").readline())
except IOError:
    open("bestfit117","w").write(str([1.0]*117))
    print "File bestfit117 not found, creating one from scratch, exitting."
    sys.exit()

print "fitting"
ndata=sum(data.flat)
def f(params):
    m=utils.normalize(iso.computeCMD(params,isos),ndata)
    return utils.loglikelihood(m,data)

d = numarray.maximum(data,1e-20)
llhC=sum( (d*numarray.log(d)).flat )
def b(params,value,iter):
    m=iso.computeCMD(params,isos)
    C=sum(data.flat)/sum(m.flat)
    x=numarray.array(params)
    file("bestfit117","w").write(str(list(x*C)))
    print "henry:",value,"tom:",2.0*(value+llhC),"iter:",iter,"norm:",C

print optimization.minc(optimization.fmin_simplex,f,pars,b)

print "done"
