""" Functions for handling isochrones (*.fits)
    2006 Henry Ferguson and Ondrej Certik, STScI.

    The isochrones are fits files, with a naming convention, e.g.:
        fm0832a110g340.fits -- [Fe/H] = -0.832, Age = 11.0 Gyr
          the 340 is a running index.

    After reading (readisos), isochrones are stored in a fairly ugly structure:
    a list of tuples, each of which has 
       (scaled metallicity, scaledage, isochrone)
    where scaled metallicity is an integer equal to [Fe/H] *1000,
    scaled age is in uints of 10**8 years, and the isochrone is an array.
"""

__version__ = "1.0.0"
__author__ = "H. C. Ferguson & O. Certik, STScI"

import glob
import os

from numpy import array
import pyfits

import utils

#SOLAR = 0.0188
SOLAR = 0.02

def readfits(filename):
    """reads the fits file 'filename', returns a 2D array of floats"""
    f = pyfits.open(filename)
    d = f[0].data
    f.close()
    #convert NumArray to NumPy
    return array(d)

def writefits(filename,data):
    """writes the 2D array 'data' to the fits file 'filename'"""
    os.remove(filename)
    pyfits.PrimaryHDU(data).writeto(filename)

def readisos(isodir="m31iso"):
    """read isochrones from 'isodir'.
    
    returns a list of (scaledfeh,scaledage,iso) for each isochrone in the order
    of 'grp' number in the filename.
    where scaled metallicity is an integer equal to [Fe/H] *1000,
    scaled age is in uints of 10**8 years, 
    and the isochrone is an array.
    """
    files = glob.glob(isodir+'/f*a*g*.fits')
    data=[0]*len(files)
    for file in files:
        f = os.path.basename(file)
        scaledfeh = int(f[2:6])
        if f[1] == 'm': scaledfeh = -scaledfeh
        scaledage = int(f[7:10])
        grp = int(f[11:14])
        assert grp>0
#        assert grp<=117
        data[grp-1]=(scaledfeh,scaledage,array(readfits(file)))
    return data

def getfilename(grp,iso):
    """ returns the filename for a given isochrone """
    feh=iso[grp][0]
    age=iso[grp][1]
    f="fp0130a095g598.fits"
    f="f"
    if feh < 0: f+="m"
    else: f+="p"
    f+="%04d"%(abs(feh))
    f+="a%04d"%(age)
    f+="g%03d.fits"%(grp+1)
    return f

def getfehagelist(data):
    """extracts fehs and ages from the data"""
    fehlist=list(set([x[0] for x in data]))
    agelist=list(set([x[1] for x in data]))
    fehlist.sort()
    agelist.sort()
    return fehlist,agelist

def searchsorted(l,num):
    for i,x in enumerate(l):
        if num<x: return i
    return len(l)

def scaled2real(scaledfeh,scaledage):
    """ Converts scaled metallicity and age to real values """
    return (scaledfeh-0.5)/1000. , (scaledage-0.5)*1.e8

def real2scaled(feh,age):
    """feh in solar units, age in years (not log)"""
    return int(1000*feh+0.5), int(age/1.e8+0.5)

metallicities=[]
ages=[]
fehagemap={}
def getn(feh,age,isos):
    """ Figure out which isochrone corresponds to the inpu feh and age 
        and return the proper index into the isochrone list.
    """
    global metallicities,ages
    if metallicities==[]:
        metallicities,ages=getfehagelist(isos)
    scaledfeh, scaledage = real2scaled(feh,age)
    feh_i = searchsorted(metallicities,scaledfeh)
    age_i = searchsorted(ages,scaledage)
    if feh_i > 0 and feh_i < len(metallicities):
       d1 = scaledfeh-metallicities[feh_i-1]
       d2 = metallicities[feh_i]-scaledfeh
       if d1 < d2:
           feh_i -=1
    if feh_i == len(metallicities):
           feh_i -=1
    if age_i > 0 and age_i < len(ages):
       d1 = scaledage-ages[age_i-1]
       d2 = ages[age_i]-scaledage
       if d1 < d2:
           age_i -=1
    if age_i == len(ages):
           age_i -=1
    f = metallicities[feh_i]
    a = ages[age_i]
    global fehagemap
    if fehagemap=={}:
        for i,x in enumerate(isos):
            fehagemap[x[0],x[1]]=i
    return fehagemap[f,a]

def getisosweights(weights,age,metallicity,isos):
    """weights, age, metallicity ..... arrays of mass,age and metallicity

    loops through all the ages - finds the iso corresponding to the age and
    metallicity and computes the isoweight.

    results is an array of iso weights
    """
    w=[0.]*len(isos)
    for m,a,feh in zip(weights,age,metallicity):
        w[getn(feh,a,isos)]+=m
    return w

def getisosweights_gauss(weights,age,metallicity,isos,sigma):
    """ same as getisoweights, but adds Gaussian scatter """ 
    w=[0.]*len(isos)
    for m,a,feh in zip(weights,age,metallicity):
        for x in utils.frange(feh-4,feh+4,0.1):
            w[getn(x,a,isos)]+=m*utils.gauss(x,feh,sigma)
    return w

def getisosweights_vgauss(weights,age,metallicity,isos,sigma,dsigmadlogt):
    """ same as getisoweights, but adds age-dependent Gaussian scatter """ 
    w=[0.]*len(isos)
    logt0 = numarray.log10(age[0])
    for m,a,feh in zip(weights,age,metallicity):
        logt = numarray.log10(a)
        sig = sigma+dsigmadlogt*(logt-logt0)
        for x in utils.frange(feh-4,feh+4,0.1):
            w[getn(x,a,isos)]+=m*utils.gauss(x,feh,sig)
    return w

def getisosweights_sgauss(weights,age,sfr,metallicity,isos,sigma,dsigmadlogt):
    """ Same as getisoweights, but adds sfr-dependent Gaussian scatter.
        **** CAUTION: not yet tested ****
    """ 
    w=[0.]*len(isos)
    logt0 = numarray.log10(age[0])
    for m,a,feh,s in zip(weights,age,metallicity,sfr):
        ss = numarray.maximum(sfr,1.e-5)
        logs = numarray.log10(ss)
        logt = numarray.log10(a)
        sig = sigma+dsigmadlogs*(logs)
        for x in utils.frange(feh-4,feh+4,0.1):
            w[getn(x,a,isos)]+=m*utils.gauss(x,feh,sig)
    return w

def sumisos_gauss_slow(mass,age,metallicity,fehlist,agelist,isos,sigma):
    """slower version"""
    summed_iso = 0.*isos[isos.keys()[0]]
    for m,a,z in zip(mass,age,metallicity):
        feh = numarray.log10(z)-numarray.log10(SOLAR)
        for x in utils.frange(feh-4,feh+4,0.1):
            summed_iso+= utils.gauss(x,feh,sigma)*m*getiso(x,a,fehlist,agelist,isos)
    return summed_iso

def sumisos_gauss(mass,age,metallicity,fehlist,agelist,isos,sigma):
    """Faster version."""
    isoweight = {}
    for k in isos.keys():
        isoweight[k] = 0.
    for m,a,z in zip(mass,age,metallicity):
        feh = numarray.log10(z/SOLAR)
        for x in utils.frange(feh-4,feh+4,0.1):
            isoweight[getindices(x,a,fehlist,agelist)]+= utils.gauss(x,feh,sigma)*m
    summed_iso = 0.*isos[isos.keys()[0]]
    for k in isoweight.keys():
        summed_iso+=isoweight[k]*isos[k]
#    p(isoweight)
    return summed_iso

def sumisos_vgauss(mass,age,metallicity,fehlist,agelist,isos,sigma,dsigmadlogt):
    """Faster version."""
    logt0 = numarray.log10(age[0])
    isoweight = {}
    for k in isos.keys():
        isoweight[k] = 0.
    for m,a,z in zip(mass,age,metallicity):
        sig = sigma+dsigmadlogt*(logt-logt0)
        feh = numarray.log10(z/SOLAR)
        for x in utils.frange(feh-4,feh+4,0.1):
            isoweight[getindices(x,a,fehlist,agelist)]+= utils.gauss(x,feh,sig)*m
    summed_iso = 0.*isos[isos.keys()[0]]
    for k in isoweight.keys():
        summed_iso+=isoweight[k]*isos[k]
#    p(isoweight)
    return summed_iso
def computefehage(pars,isos):
    """computes the array of feh vs age for the isochrones weights: pars"""
    fehlist,agelist=getfehagelist(isos)
    #s = numarray.zeros((13,9))
    s = numarray.zeros((len(fehlist),len(agelist)))
    for i,k in enumerate(isos):
        x = searchsorted(fehlist,k[0])-1
        y = searchsorted(agelist,k[1])-1
        s[x,y]=pars[i]
    return s

def computeCMD(pars,isos):
    """computes CMD using isochrones "iso" and corresponding weights "pars"."""
    #pars = list(pars)
    #print type(pars[0])
    #stop
    #t = time.time()
    #print type(pars)
    s = 0.*isos[0][2]
    assert len(isos)==len(pars)
    for i,isochrone in enumerate(isos):
        #print type(isochrone[2]), type(pars[i])
        #s+=float(pars[i])*isochrone[2]
        s+=pars[i]*isochrone[2]
    return s
