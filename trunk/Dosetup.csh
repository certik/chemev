#!/usr/bin/env csh
cd  utils
if ($?PYTHONPATH) then
    setenv PYTHONPATH `pwd`:{$PYTHONPATH}
else
    setenv PYTHONPATH `pwd`
endif
cd optimization
setenv PYTHONPATH `pwd`:{$PYTHONPATH}
cd ../../../pymc/PyMC
setenv PYTHONPATH `pwd`:{$PYTHONPATH}
setenv PYTHONPATH  build/lib.linux-x86_64-2.4/PyMC/flib.so:{$PYTHONPATH}
cd ../../../chemev
