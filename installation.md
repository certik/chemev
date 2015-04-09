This projects depends on pymc (optional), pyfits, matplotlib (optional) and numpy (numarray). It assumes that pyfits, matplotlib and numpy (numarray) are in the standard locations. See below how to install pymc, if you need it.

Create some directory, and in this dir do all the installation steps below. The "wajig" commands are the Debian way how to install the packages - in other distros, you just need to install those packages in your way.

# Installation #

```
wajig install python-numpy python-pyfits
svn checkout https://chemev.googlecode.com/svn/trunk/ chemev --username ondrej.certik
cd chemev
./unpack
```

(The ./unpack script will unpack all the isochrones)

### Install PyMC (optional) ###

Only install PyMC if you are going to use mcmc. Otherwise just skip this step.

Go to the directory where chemev is installed and:

```
wajig install python-numpy-dev python-matplotlib python-tk
svn checkout http://pymc.googlecode.com/svn/trunk/ pymc
cd pymc
python setup.py build
```

Normally, it compiles fine, but on the goods cluster, there are some problems with the fortran and python installation and it gives this error
```
Cannot open file flib.f
gcc: flib.c: No such file or directory
```
so now you need to issue these commands manually:

_Works for Ondrej_:
```
f77 -g -Wall -fno-second-underscore -fPIC -O3 -funroll-loops -march=opteron -mmmx -m3dnow -msse2 -msse -Ibuild/src.linux-x86_64-2.4 -I/usr/stsci/pyssg/Python-2.4.3/lib/python2.4/site-packages/numpy/core/include -I/usr/stsci/pyssg/Python-2.4.3/include/python2.4 -c -c PyMC/flib.f -o build/temp.linux-x86_64-2.4/PyMC/flib.o
f77 -g -Wall -fno-second-underscore -fPIC -O3 -funroll-loops -march=opteron -mmmx -m3dnow -msse2 -msse -Ibuild/src.linux-x86_64-2.4 -I/usr/stsci/pyssg/Python-2.4.3/lib/python2.4/site-packages/numpy/core/include -I/usr/stsci/pyssg/Python-2.4.3/include/python2.4 -c -c build/src.linux-x86_64-2.4/PyMC/flib-f2pywrappers.f -o build/temp.linux-x86_64-2.4/PyMC/flib-f2pywrappers.o
cd build/temp.linux-x86_64-2.4/build/src.linux-x86_64-2.4 
g77 -g -Wall -shared PyMC/flibmodule.o fortranobject.o ../../PyMC/flib.o ../../PyMC/flib-f2pywrappers.o -lg2c -o ../../../lib.linux-x86_64-2.4/PyMC/flib.so
```

_Works for Harry_:
```
f77 -g -Wall -fno-second-underscore -fPIC -O3 -funroll-loops -march=opteron -mmmx -m3dnow -msse2 -msse -Ibuild/src.linux-x86_64-2.4 -I/usr/stsci/pyssg/Python-2.4.3/lib/python2.4/site-packages/numpy/core/include -I/usr/stsci/pyssg/Python-2.4.3/include/python2.4 -c -c PyMC/flib.f -o build/temp.linux-x86_64-2.4/PyMC/flib.o
f77 -g -Wall -fno-second-underscore -fPIC -O3 -funroll-loops -march=opteron -mmmx -m3dnow -msse2 -msse -Ibuild/src.linux-x86_64-2.4 -I/usr/stsci/pyssg/Python-2.4.3/lib/python2.4/site-packages/numpy/core/include -I/usr/stsci/pyssg/Python-2.4.3/include/python2.4 -c -c build/src/PyMC/flib-f2pywrappers.f -o build/temp.linux-x86_64-2.4/PyMC/flib-f2pywrappers.o
cd build/temp.linux-x86_64-2.4/build/src
g77 -g -Wall -shared PyMC/flibmodule.o fortranobject.o ../../PyMC/flib.o ../../PyMC/flib-f2pywrappers.o -lg2c -o ../../../lib.linux-x86_64-2.4/PyMC/flib.so
```

### Set your PYTHONPATH ###

PYTHONPATH needs to be set correctly. Execute:
```
eval `~/chemev/env.py`
```
Change the path to env.py according to your configuration. You can also put this
line to your .bashrc file.

Note: the env.py returns shell commands, that add the relevant paths to your PYTHONPATH (ensures that paths are not duplicated etc.). It automatically finds the chemev and pymc package (it assumes the pymc is installed in the same parent directory as chemev). This is necessary, because the exact path to PyMC depends on the architecture.

You can check that env.py found everything by:
```
echo $PYTHONPATH
```
it should print something like:
```
/home/ondra/chemev:/home/ondra/pymc/build/lib.linux-i686-2.4
```
If you already had some paths in PYTHONPATH, they will of course be appended to it.
The first path is to the chemev package, the second is to the pymc package. If the second path is missing, it means the env.py didn't find PyMC.

Deprecated way:
```
cd chemev
source Dosetup
```

# Tests: example programs #

Now you should have directories "pymc" (optional) and "chemev" in the installation directory.

### Test simplex ###

This will run 117 isochrones fitting using the simplex method.

```
cd calc/halosimplex/
python 117isofit.py
```

This creates a file bestfit117 initialized with ones.

```
python 117isofit.py
```

This will start the fitting procedure. It will start with values in bestfit117 and everytime it finds a better solution, it will write it to bestfit117, which means that you can kill it anytime with CTRL-C and start it again.

### Test simple models ###

This will run 696 isochrones simple model fitting using the simplex method.

```
cd calc/halo_gaussian           #or "cd calc/halo_nongaussian", its just a different simple model
python fit.py
```

This will start the fitting procedure. It will start with values in bestfit and everytime it finds a better solution, it will write it to bestfit, which means that you can kill it anytime with CTRL-C and start it again.

If you want to create an initial bestfit file, then run
```
python fit.py n
```


### Test mcmc ###

You need the PyMC package to have installed (see above).
This will run 117 isochrones fitting using the mcmc method:

```
cd chemev/calc/halomcmc
python mcmc.py > /tmp/d.txt
```

after a while you can kill it with CTRL-C and inspect the file /tmp/d.txt. You should see something like this:

```
reading isochrones
fitting

Iteration 0 at 0.000184059143066

Marginal Posterior Statistics
==================================================
Node: AIC

95% HPD interval = [234.0, 234.0]
mc error = 0.0
mean = 234.0
n = 5000
quantiles =
    2.5:    234.0
    25: 234.0
    50: 234.0
    75: 234.0
    97.5:   234.0
standard deviation = 0.0

==================================================

...
```