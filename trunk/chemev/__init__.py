import os

def getisodir():
    root = os.path.dirname(__path__[0])
    return os.path.join(root,"isochrones")

isodir = getisodir()
