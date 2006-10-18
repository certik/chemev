import amoeba
import DESolver
import scipyoptimize

def f(var): 
    x=var[0]
    v=1./14*(x+4)*(x+1)*(x-1)*(x-3)+0.5
    return v

#print DESolver.fmin_de(f,[-10],[0])
#print amoeba.fmin_simplex(f,[-10],[0])
