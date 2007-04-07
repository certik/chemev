import numpy
from math import sqrt
import time

count = 0
def approx_fprime(xk,f):
    """
    Calculates a gradient of f(xk) numerically.

    """
    epsilon = sqrt(numpy.finfo(float).eps)
    t = time.time()
    #print "--1",len(xk)
    f0 = apply(f,(xk,))
    grad = numpy.zeros((len(xk),), float)
    ei = numpy.zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = epsilon
        grad[k] = (apply(f,(xk+ei,)) - f0)/epsilon
        ei[k] = 0.0
    t = time.time() - t 
    global count
    count +=1
    print "--grad %d:"%count,t
    #print xk
    return grad

def fmin_bfgs(f, x0, iter=100,callback=None):
    """Calls scipy L-BFGS optimizer (unconstrained).
    
    This is a superior unconstrained optimizer. Use it with our functions
    minmax or minc constrains to get constrained optimizing.

    TODO: there is one unnecessary call to the f(xk) - one in mycallback() and
    second in myprime() - so joining these two methods would spare us one call,
    i.e. reducing 119 calls to 118 calls per BFGS iteration.
    """
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not available"

    def mycallback(xk):
        if callback:
            callback(xk, f(xk), -1)

    def myprime(xk):
        return approx_fprime(xk, f)

    #optimize.fmin_bfgs(f, x0, fprime=myprime,
    #        maxiter = iter, callback = mycallback)
    #scipy's implementation of the bfgs converges much slower than the fortran's
    x,f,g = optimize.fmin_l_bfgs_b(f, x0, fprime = myprime,
            iprint=1, maxfun=iter, factr=10.0)
    return x

def fmin_scipy_l_bfgs_b(f, x0, iter=100):
    """Calls scipy L-BFGS-B optimizer (constrained).
    
    This function works, but for some reason, using our own constrains using the
    logistics function converges faster (in my tests by a factor of 7), so
    this function is here just for testing purposes. Use fmin_bfgs()
    instead.
    """
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not be available"

    def mycallback(xk):
        if callback:
            callback(xk, f(xk), -1)

    def myprime(xk):
        return approx_fprime(xk, f)

    b = [(0,None)]*117
    #b = [(0,1e7)]*117
    x,f,g = optimize.fmin_l_bfgs_b(f, x0, fprime = myprime, bounds = b, #approx_grad=True,
            iprint=1, maxfun=iter)
    return x,f,g

def fmin_anneal(f, x0, iter=100,callback=None):
    """Calls scipy simulated annealing optimizer (unconstrained).
    
    """
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not available"

    global bestsol, it
    bestsol = (x0, f(x0))
    it = 0
    def myf(x):
        global it, bestsol
        it += 1
        r= f(x)
        if r < bestsol[1]:
            bestsol = (x, r)
            print "%d:"%it,r
        return r

    M=1.e7
    M = 5
    min=[-M]*len(x0)
    max=[M]*len(x0)

    x, retval = optimize.anneal(myf, x0, 
            lower = min, upper = max, maxiter = iter, schedule = "fast")
    print "retval", retval
    return x
