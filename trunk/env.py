#! /usr/bin/env python

#prints shell commands that sets PYTHONPATH to correct value.
#usage: from any directory execute:
#  eval `/path/to/chemev/env.py`
#you can add this line to your .bash script

import sys, os, os.path

progpath = sys.argv[0]
#packagedir .... the absolute path to the directory, where env.py is located
#thus the path to env.py is os.path.join(packagedir,"env.py")
packagedir = os.path.abspath(os.path.dirname(progpath))

rootdir = packagedir

def prepend_path(name, value):
    sep = os.path.pathsep
    curpath = os.environ.get(name, '')
    newpath = [value] + [ x for x in curpath.split(sep) if x and x != value ]
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

print prepend_path('PYTHONPATH', rootdir)
