def fmin_scipy_bfgs(f, x0, iter=100,callback=None):
    """Calls scipy BFGS optimizer."""
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not be available"

    def mycallback(xk):
        if callback:
            callback(xk, f(xk), -1)

    optimize.fmin_bfgs(f, x0, maxiter = iter, callback = mycallback)

def fmin_scipy_l_bfgs_b(f, x0, iter=100,callback=None):
    """Calls scipy L-BFGS-B optimizer."""
    try:
        from scipy import optimize
    except ImportError:
        raise "SciPy is not installed, SciPy optimizers are not be available"

    def mycallback(xk):
        if callback:
            callback(xk, f(xk), -1)

    optimize.fmin_l_bfgs_b(f, x0, approx_grad=True, iprint=1, callback = mycallback)
