# This program assumes that pymc and flib.so are in PYTHONPATH
# Find both in the PyMC directory tree
import glob
#path= glob.glob("../../../*/*/*/*/flib.so")[0]
#path=path[:path.rfind("/")]
#path=path[:path.rfind("/")]
import sys
#sys.path.append(path)
import PyMC

def fmin_mcmc(f,x0,callback=None):
    class MySampler(PyMC.MetropolisHastings):
        def __init__(self,param_init_values,callback):
            #        print len(param_init_values),param_init_values
            PyMC.MetropolisHastings.__init__(self)
            for i,par in enumerate(param_init_values):
                self.parameter(name='p%d'%i, init_val=par)
            self.npar=len(param_init_values)
            self.c=callback
        def calculate_likelihood(self):
            params=[]
            for i in range(self.npar):
                eval("params.append(float(self.p%d))"%(i))
            return -f(params)
        def callback(self,pars,value,iter):
            params=[]
            for i in range(self.npar):
                eval("params.append(float(self.p%d))"%(i))
            self.c(params,value,iter)
    sampler = MySampler(x0,callback)
    sampler.sample(iterations=10000, burn=5000, thin=1, plot=False)
