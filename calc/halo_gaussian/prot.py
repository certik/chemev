#! /usr/bin/env python
"""Writes all the information about the best fit to the file "protocol" """

import sys
sys.path.append("/home/ondrej/py")

import utils
import iso
from params import parameters
import fit

gauss=True

isodir="/home/ondrej/data/isochrones/696/halo"
data=iso.readfits(isodir+"/datarr.fits")
isos = iso.readisos(isodir)
t=utils.frange(8,10.25,0.001)
p=parameters()
p.load()
m=fit.metallicity(t,p)
s=fit.sfr(t,p)
w=utils.calculateweights(t,s)
if gauss:
    isow=iso.getisosweights_gauss(w,10.**t,m,isos,p.sigma)
else:
    isow=iso.getisosweights(w,10.**t,m,isos)
model=iso.computeCMD(isow,isos)
#model=isos[0][2]*0.0+1.0
model=utils.normalize(model,sum(data.flat))
#model=model/sum(model.flat)

def plot_residuals(d,m):
    import pylab
    pylab.subplot(131)
    pylab.imshow(d,origin='lower',interpolation="nearest")
    pylab.subplot(132)
    pylab.imshow(m,origin='lower',interpolation="nearest")
    pylab.subplot(133)
    #pylab.imshow((d-m)/m,origin='lower',interpolation="nearest")
    pylab.imshow(d-m,origin='lower',interpolation="nearest")
    #    pylab.cool()
    pylab.savefig("graph.eps")
#    pylab.show()

plot_residuals(data,model)


f=open("protocol","w")
f.write('best fit for "%s".\n'%isodir)
if gauss:
    f.write("simple model, gaussian spreading\n")
else:
    f.write("simple model, non gauss\n")
f.write("likelihood: %r, toms: %r\n"%(utils.loglikelihood(model,data),
    utils.tomslikelihood(model,data)))
f.write("parameters:\n\n")
p.write(f)
f.write("-"*60+"\n\n")

f.write("age, [Fe/H], sfr, weight\n\n")
for x in zip(list(t),list(m),list(s),list(w)):
    f.write("%f %f %f %.10f\n"%x)
f.write("-"*60+"\n\n")

f.write("isochrone_file_name, weight\n\n")
for i,x in enumerate(isow):
    f.write("%s  %e\n"%(iso.getfilename(i,isos),x))
f.write("-"*60+"\n\n")
