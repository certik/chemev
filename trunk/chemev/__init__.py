import os

def getisodir():
    """
    Returns a path to the dir with isochrones.
    """
    root = os.path.dirname(__path__[0])
    return os.path.join(root,"isochrones")

isodir = getisodir()
