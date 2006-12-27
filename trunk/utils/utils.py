#utility functions for the project, which can be reused elsewhere.

import numarray
import math

def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."
    if end == None:
        end = start + 0.0
        start = 0.0
    if inc == None:
        inc = 1.0
    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end: break
        elif inc < 0 and next <= end: break
        L.append(next)
    return numarray.array(L)

def normalize(model,ndata):
    """ Return a model where the integral in the unmasked region is
        equal to ndata

        Arguments:
        model - masked array of unnormalized model values 
    """
    total = sum(model.flat)
    return model*ndata/total

def loglikelihood(model,data):
    """Return the -log(likelihood), i.e. the lower, the better

       Arguments:
       model ... 2D masked model array (pdf normalized to same total number
                    as for data)

       data  ... 2D masked data array (integers -- number of points per bin)

       the probability of drawing the data from the given model is:
       P=model[1,1]**data[1,1] * model[1,2]**data[1,2] * model[1,3]**data[1,3] *
       ... multiply over all pairs (i,j).

       the P is also the likelihood of the model for the given measured data.
       the higher P, the better.
       For numerical reasons, we use -log(P):
       -log(P)=-sum over (i,j) of data[i,j]*log(model[i,j])
    """
    floor=1.e-50; model = numarray.maximum(model,floor)
    l = data*numarray.log(model) # l_ij = data_ij * log(model_ij)
    return -sum(l.flat) #sums all the elements of the 2D array "l"

def tomslikelihood(model,data):
    """l=2.*(model-data+data*numarray.log(data/model))
    return sum(l.flat)
    """
    data = numarray.maximum(data,1e-20)
    C=data*numarray.log(data)
    return 2.0*(loglikelihood(model,data)+sum(C.flat))

def tomslikelihood_orig(model,data):
    data = numarray.maximum(data,1e-20)
    model = numarray.maximum(model,1e-20)
    l=2.*(model-data+data*numarray.log(data/model))
    return sum(l.flat)

def gauss(x,m,sigma):
    sigma=float(sigma)
    return math.exp(-((x-m)/sigma)**2/2.)/(sigma*math.sqrt(2*math.pi))

def getdt(logt):
    """ Return an array of time spacings of bins centered on t,
        running to the midpoints between successive t values. 

        the input is a log t. output is a real time interval (not a log).
    """
    x = logt #example: logt=[1,2,3,4]
    n = len(x)
    slope = x[1]-x[0]
    x0 = -1.*slope+x[0]  #x0=0
    slope = x[n-1]-x[n-2]
    xn = slope + x[n-1]  #xn=5
    xlo = numarray.concatenate((numarray.array([x0]),x)) #xlo=[0,1,2,3,4]
    xhi = numarray.concatenate((x,numarray.array([xn]))) #xhi=[1,2,3,4,5]
    midpts = (xhi+xlo)/2.  #midpts=[ 0.5,  1.5,  2.5,  3.5,  4.5]
    linear_midpts = 10.**midpts #=[3.16, 31.62, 316.22, 3162.27, 31622.77]
    dx = abs(linear_midpts[1:]-linear_midpts[:-1]) #[28.4,284.6,2846.0,28460.4]
    return dx

def calculateweights(t,sfr):
    #sfr is the density - SFR*dt is the total amount of new stars in the 
    #interval dt
    #get the bins centered at t - t is logage, dt contains the
    #time interval (in age, not logage) around each t
    dt = getdt(t) 
    #sfr * dt = total amount of stars in the time interval (bin) around each age
    w0=sfr*dt   
    #normalize - sum(w) = 1, so the total amount of stars in the field is 1.0
    #w(t) is the amount of new stars in the time bin around t, so
    #the sum w(t) over all t is the total amount of stars. 
    w = w0/sum(w0) 

#    sfrate = w/dt
#    sfrate = sfrate/sfrate.mean()
    return w

def plot(data2):
    import pylab
    pylab.imshow(data2,origin='lower',interpolation="nearest")
#    pylab.xticks(pylab.arange(40), [str(n+10) for n in pylab.arange(40)])
    pylab.show()

def plot2(data1,data2):
    import pylab
    pylab.subplot(121)
    pylab.imshow(data1,origin='lower',interpolation="nearest")
    pylab.subplot(122)
    pylab.imshow(data2,origin='lower',interpolation="nearest")
#    pylab.xticks(pylab.arange(40), [str(n+10) for n in pylab.arange(40)])
    pylab.show()

def plotfehsfr(t,feh,sfr):
    import pylab
    pylab.subplot(121)
    pylab.plot(t,feh)
    pylab.subplot(122)
    pylab.plot(t,sfr)
    pylab.show()

def plotfehsfr2(t,feh,sfr1,sfr2):
    import pylab
    pylab.subplot(131)
    pylab.plot(t,feh)
    pylab.subplot(132)
    pylab.plot(t,sfr1)
    pylab.subplot(133)
    pylab.plot(t,sfr2)
    pylab.show()

def p2c(r,phi):
    return r*math.cos(phi),r*math.sin(phi)