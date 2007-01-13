#! /usr/bin/env python
"""Makes summary plots of all the information the best fit. Name
   starts with optional command-line argument 'arg':
   arg_fig1 - CMDs & residuals
   arg_fig1a -- CMD
   arg_fig1b -- Model CMD
   arg_fig1c -- Residual CMD
   arg_fig2 -- sfr & metallicity dots
   arg_fig3 -- sfr & metallicity contours
   arg_fig3b -- sfr & metallicity image from isochrone weights
   arg_fig3c -- Obsolete: sfr & metallicity image from age weights
   arg_fig4 -- sfr vs time
   arg_fig5 -- metallicity vs time
   arg_fig6 -- differential metallicity distribution
   arg_fig6b -- cumulative metallicity distribution
   arg_fig7 -- differential log(age) distribution
   arg_fig7b -- differential age distribution
   arg_fig7c -- cumulative age distribution
"""

__version__ = "1.0.0"
__author__ = "H. C. Ferguson & O. Certik, STScI"

import random
import sys
sys.path.append(".")
import utils
from utils import p2c
import iso
from params import parameters
import fit
import numarray
import MLab
import LinearAlgebra
import pylab
import griddata
import Numeric

isodir="/devel/goods16/ferguson/chemev/halo.hires.b10"
#isodir="/home/ondrej/data/isochrones/696/halo"


def readin():
    """ read in the best fit:
        p,t,s,m,w,model,data,isos,isow = readin()
        parameters, time, sfr, metallicity, weights, model, data, ...
    """
    log_tmax = Numeric.log10(13.7e9)
    data=iso.readfits(isodir+"/datarr.fits")
    isos = iso.readisos(isodir)
    t=utils.frange(8,log_tmax,0.001)
    p=parameters()
    p.load()
    m,s,w,isow,model = compute_model(t,p,isos,data)
    return p,t,s,m,w,model,data,isos,isow

def compute_model(t,p,isos,data):
    """ Returns m,s,w,isow,model """
    m=fit.metallicity(t,p)
    s=fit.sfr(t,p)
    w=utils.calculateweights(t,s)
    if not p.pars.has_key('dsigmadlogt'):
        p.set('dsigmadlogt',0,0,False)
    if p.sigma > 0.:
        if p.dsigmadlogt == 0.:
            isow=iso.getisosweights_gauss(w,10.**t,m,isos,p.sigma)
        if p.dsigmadlogt != 0.:
            print "Gaussian sigma, ds/dlogt ",p.sigma,p.dsigmadlogt
            isow=iso.getisosweights_vgauss(w,10.**t,m,isos,p.sigma,p.dsigmadlogt)
    else:
        isow=iso.getisosweights(w,10.**t,m,isos)
    model=iso.computeCMD(isow,isos)
    model=utils.normalize(model,sum(data.flat))
    d = numarray.maximum(data,1e-20)
    llhC=sum( (d*numarray.log(d)).flat )
    value=utils.loglikelihood(model,data)
    print "henry:",value,"tom:",2.0*(value+llhC)
    return m,s,w,isow,model

def isovalues(isos,isow):
    x = numarray.zeros(len(isos))*0.  # age
    y = numarray.zeros(len(isos))*0.  # metallicity
    z = numarray.zeros(len(isos))*0.  # metallicity
    for i in range(len(isos)):
        feh,age = iso.scaled2real(isos[i][0],isos[i][1])
        y[i] = feh
        x[i] = numarray.log10(age)
        z[i] = isow[i] 
    return x,y,z

def fig2(isos,isow):
    """ Plot Isochrone weights, with varying point sizes """
    x,y,z = isovalues(isos,isow)
    zz = z/sum(z.flat)
    maxzz = max(zz) 
    zznorm = 1000*zz/maxzz
    pylab.scatter(10.**x/1.e9,y,zznorm,zznorm,faceted=False,alpha=0.9)
    pylab.xlim(0.5,15)
    pylab.ylim(-2.4,0.6)
    pylab.xlabel('Age (Gyr)')
    pylab.ylabel('[Fe/H]')
    return x,y,z

def fig3(isos,isow):
    """ Plot contour diagram of isochrone weights """
    x,y,z = isovalues(isos,isow)
    zz = z/sum(z.flat)
#   xi,yi = pylab.meshgrid(pylab.linspace(8.8,10.3,100),pylab.linspace(-2.5,0.9,100))
#   zi = griddata.griddata(x,y,zz,xi,yi)
    xi,yi = pylab.meshgrid(pylab.linspace(0.5,15,100),pylab.linspace(-2.5,0.9,100))
    zi = griddata.griddata(10.**x/1.e9,y,zz,xi,yi)
    pylab.contourf(xi,yi,zi,20)
    pylab.xlabel('Age (Gyr)')
    pylab.ylabel('[Fe/H]')
    return x,y,z


def fig3b(t,m,w,p):
    """ Plot image of age,metallicity weights -- nearest neighbor interpolation"""
    NREP = 50
    xi = pylab.linspace(0.1,15,512)
    yi = pylab.linspace(-2.5,0.9,512)
    z = numarray.zeros((len(xi),len(yi)),numarray.Float) # create t, metallicity array
    y = numarray.repeat(m,NREP)
    dt = t[1]-t[0]
    x = numarray.arange(t[0],t[-1]+2*dt,dt/NREP)
    x = x[0:len(y)]
    x = 10.**x/1.e9
    weight = numarray.repeat(w,NREP)
    # Find the indices in the array
    xindex = numarray.searchsorted(xi,x)
    print "shape(x), shape(y), shape(weight)", numarray.shape(x),numarray.shape(y), numarray.shape(weight)
    if p.sigma > 0:
        if p.dsigmadlogt == 0.:
            for i in range(len(y)):
                nstars = weight[i]*normgauss(yi,y[i],p.sigma)
                j = xindex[i]
                z[j,:] += nstars
        if p.dsigmadlogt != 0.:
            for i in range(len(y)):
                logt0 = numarray.log10(x[0])
                logt = numarray.log10(x[i])
                sigma = p.sigma + p.dsigmadlogt*(logt-logt0)
                nstars = weight[i]*normgauss(yi,y[i],sigma)
                j = xindex[i]
                z[j,:] += nstars
    else:
        sigma = 0.01
        for i in range(len(y)):
            nstars = weight[i]*normgauss(yi,y[i],sigma)
            j = xindex[i]
            z[j,:] += nstars
#       yindex = numarray.searchsorted(yi,y)
#       z[xindex,yindex] = z[xindex,yindex] + weight # increment the 2-d array
    zz = numarray.transpose(z)
    pylab.imshow(zz,extent=[0.1,15,-2.5,0.9],aspect='auto')
    pylab.xlabel("Age (Gyr)")
    pylab.ylabel("[Fe/H]")
    return xi,yi,zz

def fig3e(t,m,w,p):
    """ Plot image of log(age),metallicity weights -- nearest neighbor interpolation"""
    xi = pylab.linspace(8,10.25,512)
    yi = pylab.linspace(-2.5,0.9,512)
    z = numarray.zeros((len(xi),len(yi)),numarray.Float) # create t, metallicity array
    x = t
    y = m
    # Find the indices in the array
    xindex = numarray.searchsorted(xi,x)
    if p.sigma > 0:
        for i in range(len(y)):
            nstars = w[i]*normgauss(yi,y[i],p.sigma)
            print "shape(z),len(nstars)", numarray.shape(z), len(nstars)
            z[xindex[i],:] += nstars
    else:
        yindex = numarray.searchsorted(yi,y)
        z[xindex,yindex] = z[xindex,yindex] + weight # increment the 2-d array
    zz = numarray.transpose(z)
    pylab.imshow(zz,extent=[8,10.25,-2.5,0.9],aspect='auto')

def fig4(t,s): # SFR vs. time
    """ Plot SFR vs age """
    pylab.plot(10.**t/1.e9,s)
    pylab.xlim(0.1,15)
#   pylab.ylim(0.0,1.1)
    pylab.xlabel("Age (Gyr)")
    pylab.ylabel("Star-formation Rate") 

def fig4b(t,s,p):
    """ Plot SFR vs log(age) """
    pylab.plot(t,s)
    pylab.xlabel("log Age (Gyr)")
    pylab.ylabel("Star-formation Rate") 
    sfr_bezier_points(p)
    pylab.ylim(-0.1,1.5)

def fig5(t,m): # Metallicity vs. time
    """ Plot SFR vs age """
    tt = numarray.compress(t>p.s0x,t)
    mm = numarray.compress(t>p.s0x,m)
    tt = numarray.compress(tt<p.s2x,tt)
    mm = numarray.compress(tt<p.s2x,mm)
    pylab.plot(10.**tt/1.e9,mm)
    pylab.xlim(0.1,15)
    pylab.xlabel("Age (Gyr)")
    pylab.ylabel("Metallicity") 

def fig5b(t,m,p): # Metallicity vs. time
    """ Plot SFR vs age """
    tt = numarray.compress(t>p.s0x,t)
    mm = numarray.compress(t>p.s0x,m)
    tt = numarray.compress(tt<p.s2x,tt)
    mm = numarray.compress(tt<p.s2x,mm)
    pylab.plot(tt,mm)
    metallicity_bezier(p)
    pylab.xlabel("Age (Gyr)")
    pylab.ylabel("Metallicity") 

def fig6(m,w): # Metallicity distribution
    """ fig6(t,w) -- differential metallicity distribution """
    mmin,mmax,mstep = -2.5,1.0,0.2
    mmin = mmin-mstep
    mbins = numarray.arange(mmin,mmax,mstep)
    indices = numarray.searchsorted(mbins,m)
    weight = numarray.zeros(len(mbins))*0.
    for i in range(len(m)):
        weight[indices[i]] += w[i] 
    pylab.bar(mbins[1:],weight[:-1],width=mstep,edgecolor=None)
    pylab.xlabel("[Fe/H]")
    pylab.ylabel("N") 
    return mbins[1:],weight[:-1]

def fig6b(m,w): # Cumulative Metallicity distribution
    """ Plot differential metallicity distribution """
    wnorm = w/sum(w)
    i = numarray.argsort(m)
    mm = m[i] 
    ww = wnorm[i]
    pylab.plot(mm,numarray.cumsum(ww))
    pylab.xlabel("[Fe/H]")
    pylab.ylabel("N < [Fe/H]")

def fig6c(xi,yi,zz):
    """ Plot cumulative metallicity distribution """
    feh = numarray.sum(zz,1)
    pylab.plot(yi,feh)
    pylab.xlabel('[Fe/H]')
    pylab.ylabel('N')

def fig7(t,w): # Age distribution
    """ fig7b(t,w) -- differential logage distribution """
    pylab.plot(t,w)
    pylab.xlabel('Log(Age)')
    pylab.ylabel('N( Log(Age))')

def fig7b(t,w): # Age distribution
    """ fig7b(t,w) -- differential age distribution """
    tmin,tmax,tstep = 0.0,14.5,0.5
    tmin = tmin-tstep
    tbins = numarray.arange(tmin,tmax,tstep)
    t_gyr=10.**t/1.e9
    indices = numarray.searchsorted(tbins,t_gyr)
    weight = numarray.zeros(len(tbins))*0.
    for i in range(len(w)):
        weight[indices[i]] += w[i] 
    pylab.bar(tbins[1:],weight[:-1],width=tstep,edgecolor=None)
    pylab.xlabel("Age (Gyr)")
    pylab.ylabel("N") 
    return tbins[1:],weight[:-1]

def fig7c(t,w): # Cumulative Age distributions
    """ fig7c(t,w) -- cumulative age distribution """
    wnorm = w/sum(w)
    t_gyr = 10.**t/1.e9
    pylab.plot(t_gyr,numarray.cumsum(wnorm))
    pylab.xlabel("Age (Gyr)")
    pylab.ylabel("N < Age")

def sfr_bezier_points(p):
    """ Plot the Bezier points for sfr(t) """
    s0tx = p.s0x + p.s0tx*(p.s1x-p.s0x)
    s0ty = p.s0y + p.s0ty*(p.s1y-p.s0y)
    s1tx = p.s1x - p.s1tx*(p.s2x-p.s1x)
    s1ty = p.s1y + p.s1ty*(p.s1y-p.s0y)
    s1btx = p.s1x + p.s1tx*(p.s2x-p.s1x)
    s1bty = p.s1y - p.s1ty*(p.s1y-p.s0y)
    s2tx = p.s2x - p.s2tx*(p.s2x-p.s1x)
    s2ty = p.s2y + p.s2ty*(p.s1y-p.s2y)

    p0 = p.s0x,p.s0y,"p0"
    t0 = s0tx,s0ty,"t0"
    p1 = p.s1x,p.s1y,"p1"
    t1 = s1tx,s1ty,"t1"
    t1b = s1btx,s1bty,"t1b"
    p2 = p.s2x, p.s2y,"p2"
    t2 = s2tx,s2ty,"t2"

    points = zip(p0,t0,t1,p1,t1b,t2,p2)
    for i in range(len(points[0])):
        print points[0][i],points[1][i],points[2][i]
    pylab.plot(points[0],points[1])
    pylab.plot(points[0],points[1],'ro')
    points = zip(t1b,t2,p2)
    pylab.plot(points[0],points[1])
    pylab.plot(points[0],points[1],'go')

def metallicity_bezier(p):
    "Fe/H vs age with Bezier control points"
    m0y =  p.m0y*(p.m0ymax - p.m1y) + p.m1y # Metallicity increase with time
    m1tx = p.s2x-p.m1tx*(p.s2x-p.s0x)
    m1ty = p.m1y+p.m1ty*(m0y-p.m1y)
    m0tx = p.s0x+p.m0tx*(p.s2x-p.s0x)
    m0ty = m0y-p.m0ty*(m0y-p.m1y)

    p0 = p.s0x,m0y,"m0"
    t0 = m0tx,m0ty,"t0"
    p1 = p.s2x,p.m1y,"m1"
    t1 = m1tx,m1ty,"t1"

    points = zip(p0,t0,t1,p1)
    for i in range(len(points[0])):
        print points[0][i],points[1][i],points[2][i]
    pylab.plot(points[0],points[1])
    pylab.plot(points[0],points[1],'ro')

def plot_residuals(d,m):
    """ Plot data, model, and residuals """
    extent = [-0.9,0.1,30.5,26.5]
    pylab.subplot(131)
    pylab.imshow(d,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.ylabel("I (STMAG)")
    pylab.title("Data")
    pylab.subplot(132)
    pylab.imshow(m,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.yticks([])
    pylab.title("Model")
    pylab.xlabel("R-I (STMAG)")
    pylab.subplot(133)
    #pylab.imshow((d-m)/m,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.imshow(d-m,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.yticks([])
    pylab.title("Residuals")
    pylab.subplots_adjust(left=0.1,right=0.95,wspace=0.02)
#    pylab.show()

def gauss(x,m,sigma):
    """ Return a Gaussian with mean, sigma for a numarray x """
    sigma2 = sigma*sigma
    exponent = -(x-m)**2/(2.*sigma2)
    exponent = numarray.choose(exponent<-700.,(exponent,-700.))
    try:
       result = numarray.exp(exponent)/numarray.sqrt(2.*numarray.pi*sigma2)
    except OverflowError:
       print "gauss: overflow error"
       print "sigma = ", sigma
       print "m = ", m
       print "x = ", x
       sys.exit()
    return result

def normgauss(x,m,sigma):
    g = gauss(x,m,sigma)
    return g/sum(g)

def changenode(node,value,t,p,isos,data):
    """ Change a node value and replot the cmd in figure 1 and fig4b in
        figure 2 """
    p.set(node,value,p.getmin(node),p.getmax(node),True)
    m,s,w,isow,model = compute_model(t,p,isos,data)
    pylab.figure(1)
    extent = [-0.9,0.1,30.5,26.5]
    pylab.subplot(132)
    pylab.imshow(model,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.yticks([])
    pylab.title("Model")
    pylab.xlabel("R-I (STMAG)")
    pylab.subplot(133)
    #pylab.imshow((d-m)/m,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.imshow(data-model,origin='lower',interpolation="nearest",extent=extent,aspect=0.6)
    pylab.yticks([])
    pylab.title("Residuals")
    pylab.figure(2)
    pylab.clf()
    fig4b(t,s,p)
    return p,m,s,w,isow,model

if __name__ == "__main__":
    if len(sys.argv) > 1:
        root = sys.argv[1]+"_"
    else:
        root=""
    p,t,s,m,w,model,data,isos,isow = readin()
    pylab.rc('font',family='sans-serif')
    plot_residuals(data,model)
    pylab.rc('axes',linewidth=1.5)
    pylab.rc('axes',titlesize='larger')
    pylab.rc('axes',labelsize='larger')
    pylab.rc('xtick',labelsize='large')
    pylab.rc('ytick',labelsize='large')
    pylab.rc('lines',linewidth=3.0)
    pylab.savefig(root+"fig1.png")
    pylab.clf(); fig2(isos,isow)
    pylab.savefig(root+"fig2.png")
    pylab.clf(); xi,yi,zz = fig3b(t,m,w,p)
    pylab.savefig(root+"fig3b.png")
    pylab.clf(); fig4(t,s)
    pylab.savefig(root+"fig4.png")
    pylab.clf(); fig5(t,m)
    pylab.savefig(root+"fig5.png")
    pylab.clf(); fig5b(t,m,p)
    pylab.savefig(root+"fig5b.png")
#   pylab.clf(); fig6(m,w)
#   pylab.savefig(root+"fig6.png")
    pylab.clf(); fig6b(m,w)
    pylab.savefig(root+"fig6b.png")
    pylab.clf(); fig6c(xi,yi,zz)
    pylab.savefig(root+"fig6c.png")
    pylab.clf(); fig7(t,w)
    pylab.savefig(root+"fig7.png")
    pylab.clf(); fig7b(t,w)
    pylab.savefig(root+"fig7b.png")
    pylab.clf(); fig7c(t,w)
    pylab.savefig(root+"fig7c.png")
    pylab.clf(); fig4b(t,s,p)
    pylab.savefig(root+"fig4b.png")


    f=open("protocol","w")
    f.write('best fit for "%s".\n'%isodir)
    if p.sigma>0:
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
