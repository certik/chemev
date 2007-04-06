#!/usr/bin/env python

#import cgitb; cgitb.enable(format="text")
import math
from math import pi

import pylab
from numpy import array, maximum, log, sum

import chemev
from chemev import bezier, iso, utils, optimization
from chemev.utils import p2c
from chemev.params import parameters


def metallicity(t,p):
    "Fe/H vs logage on the grid 't', calculated from 'params'"
    return bezier.interpolate(bezier.curve((
            bezier.point((8,p.m0y), p2c(p.m0cr,p.m0cphi)),
            bezier.point((10.25,p.m1y),p2c(p.m1cr,p.m1cphi))
            ),8),t)

def sfr(t,p):
    "SFR vs logage on the grid 't', calculated from 'params'"
    return bezier.interpolate(bezier.curve((
            bezier.point((8,p.s0y),p2c(p.s0cr,p.s0cphi)),
            bezier.point((p.s1x,p.s1y),p2c(p.s1cr,p.s1cphi)),
            bezier.point((10.25,p.s2y),p2c(p.s2cr,p.s2cphi))
            ),8),t)

def simul(isodir):
    params=parameters()
    eps=0.01
    #feh
    params.set("m0y"   ,0.9, -1,1,True)
    params.set("m0cphi",-0.001,  -pi/2+eps,0,True)
    params.set("m0cr"  ,0.01,  0,2,True)
    params.set("m1y"   ,-0.9,  -1,1,True)
    params.set("m1cphi",1.67,  pi/2+eps,pi,True)
    params.set("m1cr"  ,1.06,  0,2,True)

#    params.set("sigma" ,0.6,  0,1, True)
    #sfr
    params.set("s0y"   ,1.11,   0.0,2,True)
    params.set("s0cphi",0.0003,   0,pi/2-eps,True)
    params.set("s0cr"  ,1,   0,1,True)
    params.set("s1x"   ,9.9,   8,10,True)
    params.set("s1y"   ,5.9, 1.5,10,True)
    params.set("s1cphi",4.61,   pi/2+eps,pi*3/2-eps,True)
    params.set("s1cr"  ,0.49,   0,0.5,True)
    params.set("s2y"   ,0.89, 0,2,True)
    params.set("s2cphi",3.04,   pi/2+eps,pi-eps,True)
    params.set("s2cr"  ,0.001,   0,2,True)
    #if len(sys.argv) == 2: #run with a param to start from the beginning
    params.save()
    params.load()

    data=iso.readfits(isodir+"/datarr.fits")
    isos = iso.readisos(isodir)
    t=utils.frange(8,10.25,0.001)
    def f(par):
        params.setvalues(par)
        w=utils.calculateweights(t,sfr(t,params))
        isow=iso.getisosweights(w,10.**t,metallicity(t,params),isos)
        #isow=iso.getisosweights_gauss(w,10.**t,metallicity(t,params),isos,
        #        params.sigma)
        m=iso.computeCMD(isow,isos)
        m=utils.normalize(m,sum(data.flat))
        return utils.loglikelihood(m,data)

    d = maximum(data,1e-20)
    llhC=sum( (d*log(d)).flat )
    def b(par,value,iter):
        params.setvalues(par)
        params.save()
        print "henry:",value,"tom:",2.0*(value+llhC),"iter:",iter

    optimization.minmax(optimization.fmin_bfgs,f,
            params.getvalues(),params.min(),params.max(),
            callback=b, iter=20)

if __name__ == "__main__":
    simul(chemev.isodir+"/696/halo")
    #pylab.show()
