import sys
#sys.path.append("/home/ondrej/data/PyMC-1.0-beta2/build/lib.linux-x86_64-2.4")
sys.path.append("/home/ondrej/data/tests/pymc/build/lib.linux-x86_64-2.4")
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
