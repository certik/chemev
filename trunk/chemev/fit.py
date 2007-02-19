#!/usr/bin/env python

#import cgitb; cgitb.enable(format="text")
""" This routine implements a simple Bezier-curve 'chemical evolution' model.
    In our current application, the star-formation rate sfr(t) is 
    constrained to have a single maximum, and the metallicity evolution
    Z(t) is constrained to increase monotonically. These restrictions can
    be relaxed without modifying the code. A description of the parameters
    is as follows:
       s0x, s0y -- time and sfr at the latest time (smallest t value)
       s1x, s1y -- time and sfr at the peak of the sfr
       s2x, s2y -- time and sfr at the earliest time
       s0tx, s0ty -- tangent control point for s0. s0tx is expressed as
          a fractional distance between s0x and s1x, s0ty is expressed 
          as a fractional distance between s0y and s1y.
       s1tx, s1ty -- tangent control point for s1, also expressed as a fraction	
       s2tx, s2ty -- tangent control point for s2, also expressed as a fraciton
       m0ymax -- maximum allowed metallicity
       m0x,m0y -- metallicity at the latest time (smallest t value); 
           m0y is expressed as a fractional distance between m1y and m0ymax
       m1x,m1y -- metallicity at the earliest time
       m0tx, m0ty -- tangent control point for m0, expresed as a fraction
       m1tx, m1ty -- tangent control point for m0, expresed as a fraction
       sigma -- Gaussian spread in metallicity
       dsigmadlogt -- change of sigma with log t
"""

__version__ = "1.0.0"
__author__ = "H. C. Ferguson & O. Certik, STScI"
    
import time
import sys
sys.path.append("../../utils")
import math
from math import pi

import pylab
import numarray

import bezier
import iso
import utils
import optimization
from utils import p2c
from params import parameters

def metallicity(t,p):
    "Fe/H vs logage on the grid 't', calculated from 'params'"
    m0y =  p.m0y*(p.m0ymax - p.m1y) + p.m1y # Metallicity increase with time
    m1tx = -p.m1tx*(p.s2x-p.s0x)
    m1ty = +p.m1ty*(m0y-p.m1y)
    m0tx = +p.m0tx*(p.s2x-p.s0x)
    m0ty = -p.m0ty*(m0y-p.m1y)
    return bezier.interpolate(bezier.curve((
            bezier.point((p.s0x,m0y), (m0tx,m0ty)),
            bezier.point((p.s2x,p.m1y),(m1tx,m1ty))
            ),8),t)

def sfr(t,p):
    "SFR vs logage on the grid 't', calculated from 'params'"
    # Set tangent points based on fractional distances between control points
    maxy02 = max(p.s2y,p.s1y)
    maxs1y = 2.*p.s1y-maxy02
    s1x = p.s0x+p.s1x*(p.s2x-p.s0x)
    s0tx = +p.s0tx*(s1x-p.s0x)
    s0ty = +p.s0ty*(p.s1y-p.s0y)
    s1ty = +p.s1ty*(maxs1y-p.s0y)
    s1tx = -p.s1tx*(p.s2x-s1x)
    s2tx = -p.s2tx*(p.s2x-s1x)
    s2ty = +p.s2ty*(p.s1y-p.s2y)
#   s1x = p.s0x+p.s1x*(p.s2x-p.s0x)

    sfr = bezier.interpolate(bezier.curve((
            bezier.point((p.s0x,p.s0y),(s0tx,s0ty)),
            bezier.point((s1x,p.s1y),(s1tx,s1ty)),
#           bezier.point((s1x,p.s1y),(s1tx,s1ty)),
            bezier.point((p.s2x,p.s2y),(s2tx,s2ty))
            ),8),t)
    sfr = numarray.where(t<p.s0x,0.,sfr)
    sfr = numarray.where(t>p.s2x,0.,sfr)
    return sfr


def simul(isodir):
    """ Read in parameters, data and isochrones. Create callback functions
    for the optimization routine, one of which will return the log(likelihood)
    and the other of which will print the best-fit parameter values. Having
    done this, call the optimization routine to minimize log(L).
    """
    log_tmax = math.log10(13.7e9)
    params=parameters()
    eps=0.01
    #feh
    params.set("m0y"   ,1.0, 0,2.5,True)
    params.set("m0cphi",-0.001,  -pi/2+eps,0,True)
    params.set("m0cr"  ,0.01,  0,2,True)
    params.set("m1y"   ,-2.5,  -2.5,0.9,True)
    params.set("m1cphi",1.67,  pi/2+eps,pi,True)
    params.set("m1cr"  ,1.06,  0,2,True)
    params.set("sigma" ,0.2,  0,1, True)
    params.set("dsigmadlogt" ,0.2,  -1,1, True)

    #sfr
    params.set("s0x"   ,8.0,   8.0,9.0,True)
    params.set("s0y"   ,0.5,   0.0,1,True)
    params.set("s0tx"  ,0.1,   0.,1.,True)
    params.set("s0ty"  ,0.1,   0,1,True)
    params.set("s1tx"  ,0.1,   0.,1.,True)
    params.set("s1ty"  ,0.1,   -1,1,True)
    params.set("s1x"   ,0.5,   0,1,True)
    params.set("s1y"   ,1.0,   0.0,1.0,False)
    params.set("s2x"   ,log_tmax,  9.5,10.25,True)
    params.set("s2y"   ,0.1,   0,1.0,True)
    params.set("s2tx"  ,0.1,   0.,1.,True)
    params.set("s2ty"  ,0.1,   0.,1.,True)

    if len(sys.argv) == 2:
        if sys.argv[1] == "start": #run with a param to start from the beginning
            params.save()
    params.load()
    if not params.pars.has_key('dsigmadlogt'):
        params.set('dsigmadlogt',0.,0,False)
    if not params.pars.has_key('dsigmadlogs'):  # Hook for SFR-depenedent spread; not fully implemented 
        params.set('dsigmadlogs',0.,0,False)
    if len(sys.argv) == 2:
        if sys.argv[1] == "nudge": #Tweak the values near their limits
             print "Nudging parameters near the limits"
             p1 = params.getl()
             utils.nudge(params)
             p2 = params.getl()
             for pp1,pp2 in zip(p1,p2):
                 if pp1[1] != pp2[1]:
                      print "%s %.8f -> %.8f" % (pp1[0],pp1[1],pp2[1])

    data=iso.readfits(isodir+"/datarr.fits")
    isos = iso.readisos(isodir)
    t=utils.frange(8,log_tmax,0.001)
    def f(par):
        params.setvalues(par)
        p = params
        w=utils.calculateweights(t,sfr(t,params))
        # isow=iso.getisosweights(w,10.**t,metallicity(t,params),isos)
        if p.sigma > 0.:
            if p.dsigmadlogt == 0.:
                isow=iso.getisosweights_gauss(w,10.**t,metallicity(t,p),isos,p.sigma)
            if p.dsigmadlogt != 0.:
#               print "Gaussian sigma, ds/dlogt ",p.sigma,p.dsigmadlogt
                isow=iso.getisosweights_vgauss(w,10.**t,metallicity(t,p),isos,p.sigma,p.dsigmadlogt)
            if p.dsigmadlogs != 0.: # Hook for SFR-depenedent spread; not fully implemented 
                isow=iso.getisosweights_sgauss(w,10.**t,sfr(t,params),metallicity(t,p),
                   isos,p.sigma,p.dsigmadlogs)
        else:
            isow=iso.getisosweights(w,10.**t,metallicity(t,p),isos)

        m=iso.computeCMD(isow,isos)
        m=utils.normalize(m,sum(data.flat))
        return utils.loglikelihood(m,data)

    d = numarray.maximum(data,1e-20)
    llhC=sum( (d*numarray.log(d)).flat )
    def b(par,value,iter):
        params.setvalues(par)
        params.save()
        print "henry:",value,"tom:",2.0*(value+llhC),"iter:",iter,time.ctime()
        sys.stdout.flush()

    optimization.minmax(optimization.fmin_simplex,f,
            params.getvalues(),params.min(),params.max(),b)

if __name__ == "__main__":
    from socket import gethostname; print gethostname()
    from time import ctime; print ctime()
    simul("/devel/goods16/ferguson/chemev/halo.hires.b10")
#   simul("/home/ondrej/data/isochrones/696/halo")
    #pylab.show()
