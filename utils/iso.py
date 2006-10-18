#functions for handling isochrones (*.fits)
#2006 Henry Ferguson and Ondrej Certik, STScI

import glob
import os

import numarray 
import pyfits

import utils

#isodir="m31iso"

#SOLAR = 0.0188
SOLAR = 0.02

def readiso(file):
    f = pyfits.open(file)
    d = f[0].data
    f.close()
    return d

def writeiso(filename,data):
    os.remove(filename)
    pyfits.PrimaryHDU(data).writeto(filename)

def readisos(isodir="m31iso"):
    files = glob.glob(isodir+'/f*a*g*.fits')
    isofile={}
    for file in files:
        f = os.path.basename(file)
        feh =  int(f[2:6])
        if f[1] == 'm': feh = -feh
        age =  int(f[7:10])
        grp =  int(f[11:14])
        isofile[feh,age] = file
    fehlist=list(set([x[0] for x in isofile]))
    agelist=list(set([x[1] for x in isofile]))
    fehlist.sort()
    agelist.sort()
    print fehlist
    print agelist
    data={}
    for f in fehlist:
        for a in agelist:
            if (f,a) not in isofile: 
                print f,a
                continue
            data[f,a] = readiso(isofile[f,a])
    return numarray.array(fehlist),numarray.array(agelist),data,isofile

def getiso(feh,age,metallicities,ages,isochrones):
    f,a=getindices(feh,age,metallicities,ages)
    return isochrones[f,a]

def searchsorted(l,num):
    for i,x in enumerate(l):
        if num<x: return i
    return len(l)

def scaled2real(scaledfeh,scaledage):
    return (scaledfeh-0.5)/1000. , (scaledage-0.5)*1.e8

def real2scaled(feh,age):
    return int(1000*feh+0.5), int(age/1.e8+0.5)

def getindices(feh,age,metallicities,ages):
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
    return f,a

def sumisos(mass,age,metallicity,fehlist,agelist,isos):
    """mass, age, metallicity ..... arrays of mass and metallicity and the 
    corresponding age.

    for each metallicity and age, there is one iso in the dir.

    loops through all the ages - finds the corresponding metallicity and thus
    iso and sums the corresponding mass * iso. 

    results is a summed iso.
    """
    summed_iso = 0.*isos[isos.keys()[0]]
    for m,a,z in zip(mass,age,metallicity):
        feh = numarray.log10(z)-numarray.log10(SOLAR)
        summed_iso += m*getiso(feh,a,fehlist,agelist,isos)
    return summed_iso

def sumisos_fast(mass,age,metallicity,fehlist,agelist,isos):
    """Faster version."""
    isoweight = {}
    for k in isos.keys():
        isoweight[k] = 0.
    for m,a,z in zip(mass,age,metallicity):
        feh = numarray.log10(z/SOLAR)
        isoweight[getindices(feh,a,fehlist,agelist)]+= m
    summed_iso = 0.*isos[isos.keys()[0]]
    for k in isoweight.keys():
        summed_iso+=isoweight[k]*isos[k]
    return summed_iso

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

def p(d):
    x=d.items()
    def cmp(x): return x[1]
    x.sort(key=cmp,reverse=True)
    for a in x[:10]:
        print scaled2real(a[0][0],a[0][1]),a[1]
