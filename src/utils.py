import numpy as np

def lowercase_array(x):
    if x.dtype.names:
        x.dtype.names = tuple(d.lower() for d in x.dtype.names)
