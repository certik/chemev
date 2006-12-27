class parameters:
    def __init__(self):
        self.pars={}
    def save(self,filename="bestfit"):
        file(filename,"w").write(str(self.pars))
    def load(self,filename="bestfit"):
        self.pars=eval(file(filename).readline())
    def set(self,par,value,min=None,max=None,fit=False):
        self.pars[par]=[value,min,max,fit]
        if min==None or max==None:
            assert min==max
            assert not fit
    def getvalues(self):
        return [self.pars[x][0] for x in self.pars.keys() if
            self.pars[x][3]]
    def min(self):
        return [self.pars[x][1] for x in self.pars.keys() if
            self.pars[x][3]]
    def max(self):
        return [self.pars[x][2] for x in self.pars.keys() if
            self.pars[x][3]]
    def getl(self):
        l=[]
        for x in self.pars: l.append(( x,self.pars[x][0]))
        def key(x): return x[0]
        l.sort(key=key)
        return l
    def write(self,f):
        l=self.getl()
        for x in l: f.write("%6s: %r\n"%(x[0],x[1]))
    def show(self):
        l=self.getl()
        for x in l: print "%6s: %r"%(x[0],x[1])
    def setvalues(self,par):
        i=0
        for x in self.pars.keys():
            if self.pars[x][3]:
                self.pars[x][0]=par[i]
                i+=1
    def __getattr__(self,name):
        if self.pars.has_key(name):
            return self.pars[name][0]
        else:
            raise AttributeError,name
