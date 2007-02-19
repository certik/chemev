# This program assumes that PyMC are in PYTHONPATH
# Try to import it. The have_mcmc variable can be later used to check if we
# succeeded.
have_mcmc=False
try:
    import PyMC
    have_mcmc=True
except ImportError, inst:
    if inst[0]=="No module named PyMC":
        #this is ok, the PyMC is not found
        pass
    else:
        #PyMC was found, but failed for some other reason, reraise
        raise

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
