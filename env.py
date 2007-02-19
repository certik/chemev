#! /usr/bin/env python

#prints shell commands that sets PYTHONPATH to correct value.
#usage: from any directory execute:
#  eval `/path/to/chemev/env.py`
#you can add this line to your .bash script

import sys, os, os.path, glob

def prepend_path(name, values):
    sep = os.path.pathsep
    curpath = os.environ.get(name, '')
    newpath = values + [ x for x in curpath.split(sep) 
            if x and (not x in values) ]
    return setenv(name, sep.join(newpath))

def setenv(name, value):
    shell = os.environ.get('SHELL', '')
    comspec = os.environ.get('COMSPEC', '')
    if shell.endswith('csh'):
        cmd = 'setenv %s "%s"' % (name, value)
    elif shell.endswith('sh'):
        cmd = '%s="%s"; export %s' % (name, value, name)
    elif comspec.endswith('cmd.exe'):
        cmd = 'set %s=%s' % (name, value)
    else:
        assert False, 'Shell not supported.'
    return cmd

def findpymc(start):
    x=os.path.join(os.path.dirname(start),"pymc")
    x = glob.glob(x+"/build/*/*/flib.so")
    if len(x) == 1: 
        return os.path.dirname(os.path.dirname(x[0]))

progpath = sys.argv[0]
#packagedir .... the absolute path to the directory, where env.py is located
#thus the path to env.py is os.path.join(packagedir,"env.py")
packagedir = os.path.abspath(os.path.dirname(progpath))
pythonpath=[packagedir]
#guess, where the pymc package is instaled
pymcdir = findpymc(packagedir)
if pymcdir:
    pythonpath+=[pymcdir]

print prepend_path('PYTHONPATH', pythonpath)
