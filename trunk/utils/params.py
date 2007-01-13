""" Implements a class for holding parameter values for chemical evolution fitting.
    A lame attempt was made to make this iterable and indexed, but it doesn't 
    entirely work (specifically not for slices). The advantage of making it indexed
    is that one could pass the entire object and not just the parameter values to
    the optimization routines. An __add__ method would take care of moving from
    one step to the next in the search of parameter space, while the callback
    routine could handle the treatment of limits: whether to penalize them
    or whatever and this would be transparent to the optimization routuine.
"""

__version__ = "1.0.0"
__author__ = "H. C. Ferguson & O. Certik, STScI"

class parameters:
    def __init__(self):
        """ Set up the class."""
        self.pars={}
        self.index = 0
    def save(self,filename="bestfit"):
        """ Save the best fit parameters to a file."""
        file(filename,"w").write(str(self.pars))
    def load(self,filename="bestfit"):
        """ Read the previous best fit from a file."""
        self.pars=eval(file(filename).readline())
    def set(self,par,value,min=None,max=None,fit=False):
        """ Set the value of a single parameter named par. """
        self.pars[par]=[value,min,max,fit]
        if min==None or max==None:
            assert min==max
            assert not fit
    def getvalues(self):
        """ Return the value of all the non-fixed parameters as a list """
        return [self.pars[x][0] for x in self.pars.keys() if
            self.pars[x][3]]
    def min(self):
        """ Return the minimum allowed value of all the non-fixed parameters 
            as a list. 
        """
        return [self.pars[x][1] for x in self.pars.keys() if
            self.pars[x][3]]
    def max(self):
        """ Return the maximum allowed value of all the non-fixed parameters 
            as a list. 
        """
        return [self.pars[x][2] for x in self.pars.keys() if
            self.pars[x][3]]
    def getl(self):
        """ Return a list of tuples: (parameter name, value).""" 
        l=[]
        for x in self.pars: l.append(( x,self.pars[x][0]))
        def key(x): return x[0]
        l.sort(key=key)
        return l
    def write(self,f):
        """ Write out current values to a file """
        l=self.getl()
        for x in l: f.write("%6s: %r\n"%(x[0],x[1]))
    def show(self):
        """ Print current values to the screen """
        l=self.getl()
        for x in l: print "%6s: %r"%(x[0],x[1])
    def setvalues(self,par):
        """ Set the values of the parameters. The input parameter par is
            a list which must be sorted in the same order as self.pars.keys().
        """
        i=0
        for x in self.pars.keys():
            if self.pars[x][3]:
                self.pars[x][0]=par[i]
                i+=1
    def getmin(self,name):
        """ Return the minimum allowed value for a single parameter. """
        if self.pars.has_key(name):
            return self.pars[name][1]
        else:
            raise AttributeError,name
    def getmax(self,name):
        """ Return the maximum allowed value for a single parameter. """
        if self.pars.has_key(name):
            return self.pars[name][2]
        else:
            raise AttributeError,name
    def __getattr__(self,name):
        """ Return the current value for a single parameter, acessed by name. """
        if self.pars.has_key(name):
            return self.pars[name][0]
        else:
            raise AttributeError,name
    def __getitem__(self,i):
        """ Return the current value for a single parameter, accessed by index. """
        return self.getl()[i][1]
    def __setitem__(self,i,value):
        """ Set the current value for a single parameter, accessed by index. """
        key = self.getl()[i][0]
        self.pars[key][0] = value
    def __iter__(self):
        return self
    def iterkeys(self):
        self.__iter__()
    def next(self):
        self.index += 1
        if self.index > len(self.pars): 
             raise StopIteration
        return self.__getitem__(self.index)
        
