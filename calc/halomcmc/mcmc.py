#import cgitb; cgitb.enable(format="text")

import sys
import pickle

import numpy as numarray

from chemev import iso,utils,optimization

if not optimization.have_mcmc:
    print "PyMC not installed correctly. Check PYTHONPATH."
    sys.exit()

print "reading isochrones"
#isodir="m31iso"
#isodir="m117halo"
isodir="../../isochrones/117/halo"
data=iso.readfits(isodir+"/datarr.fits")
isos = iso.readisos(isodir)

print "fitting"
pars=eval(file("bestfit117").readline())

ndata=sum(data.flat)
d = numarray.maximum(data,1e-20)
llhC=sum( (d*numarray.log(d)).flat )
iter=0
out=file("mcmc.out","w")
def f(params):
    global iter
    m=iso.computeCMD(params,isos)
    C=ndata/sum(m.flat)
    x=numarray.array(params)
#    out.write(str(list(x*C))+"\n")
    pickle.dump(x*C,out,protocol=-1)
    value=utils.loglikelihood(m*C,data)
#    out.write("henry: %r tom: %r iter: %r norm: %r\n" 
#        %(value,2.0*(value+llhC),iter,C))
    pickle.dump((value,2.0*(value+llhC),iter,C),out,protocol=-1)
    #out.flush()
    iter+=1
    return value

print optimization.minc(optimization.fmin_mcmc,f,pars)

print "done"
