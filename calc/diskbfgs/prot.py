#! /usr/bin/env python
"""Writes all the information about the best fit to the file "protocol" """

import sys
from chemev import iso, isodir, utils

print "reading isochrones"
data = iso.readfits(isodir+"/117/halo/datarr.fits")
isos = iso.readisos(isodir+"/117/halo")

isow = eval(open("bestfit117").readline())
model = iso.computeCMD(isow,isos)
model = utils.normalize(model,sum(data.flat))

def plot_residuals(d,m,aspect=0.2):
    """
    todo: the same color range for both data and model. the rediduals are ok
    """
    import pylab
    pylab.figure()
    pylab.subplot(221)
    pylab.imshow(d,origin='lower',interpolation="nearest",aspect=aspect)
    pylab.title("data")
    pylab.colorbar()
    pylab.subplot(222)
    pylab.imshow(m,origin='lower',interpolation="nearest",aspect=aspect)
    pylab.title("model")
    pylab.colorbar()
    pylab.subplot(223)
    residuals = d-m
    pylab.imshow(residuals,origin='lower',interpolation="nearest",aspect=aspect)
    pylab.title("data-model")
    pylab.colorbar()
    pylab.subplot(224)
    residuals = d-m
    residuals /= m
    pylab.imshow(residuals,origin='lower',interpolation="nearest",aspect=aspect)
    pylab.title("(data-model)/model")
    #pylab.cool()
    #pylab.savefig("graph.eps")
    pylab.colorbar()
    pylab.savefig("graph.png")
    #pylab.show()

def plot_weights(w,isos,aspect=0.4):
    graph=iso.computefehage(w,isos)
    import pylab
    pylab.figure()
    pylab.imshow(graph,origin='lower',interpolation="nearest",aspect=aspect)
    pylab.title("isochrones weights")
    #pylab.cool()
    #pylab.savefig("graph.eps")
    pylab.colorbar()
    pylab.savefig("graph-weights.png")
    #pylab.show()


print "likelihood: %r, toms: %r\n"%(utils.loglikelihood(model,data),
    utils.tomslikelihood(model,data))
plot_residuals(data,model)
plot_weights(isow,isos)
