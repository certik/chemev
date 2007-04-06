import numpy
from math import sqrt
import time

count = 0

def approx_fprime(xk,f):
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

def fmin_scipy_bfgs(f, x0, iter=100,callback=None):
    """Calls scipy BFGS optimizer."""
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not be available"

    def mycallback(xk):
        if callback:
            callback(xk, f(xk), -1)

    def myprime(xk):
        return approx_fprime(xk, f)

    #optimize.fmin_bfgs(f, x0, fprime=myprime,
    #        maxiter = iter, callback = mycallback)
    #scipy's implementation of the bfgs converges much slower than the fortran's
    x,f,g = optimize.fmin_l_bfgs_b(f, x0, fprime = myprime,
            iprint=1, maxfun=iter)
    return x

def fmin_scipy_l_bfgs_b_old(f, x0, iter=100,callback=None):
    """Calls scipy L-BFGS-B optimizer."""
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not be available"

    def mycallback(xk):
        if callback:
            callback(xk, f(xk), -1)

    optimize.fmin_l_bfgs_b(f, x0, approx_grad=True, iprint=1)

def fmin_scipy_l_bfgs_b(f, x0, iter=100):
    """Calls scipy L-BFGS-B optimizer."""
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
