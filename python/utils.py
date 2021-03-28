import numpy as np
from scipy.signal import convolve


def mean_with_err(arr, axis=-1):
    return arr.mean(axis=axis), arr.std(axis=axis) / np.sqrt(arr.shape[axis])


def median_with_err(arr, axis=-1):
    mu = np.median(arr, axis=axis)
    return mu, np.mean(np.abs(arr - mu), axis=axis)


def moving_average(x, w, stride=1):
    s = (x.shape[0]//stride)*stride
    x = x[-s:]
    x = x.reshape((x.shape[0]//stride, stride))
    assert w % stride == 0
    
    return convolve(x, np.ones((w//stride, stride)), 'valid', method="direct")[:,0] / w

def moving_mean_err(x, w, stride):
    mean_x_square = moving_average((x**2).astype("float64"), w, stride)
    mean_x = moving_average(x, w, stride)
    err = np.sqrt(mean_x_square - mean_x**2)/np.sqrt(w)
    return mean_x, err

def cum_mean_err(x):
    x = x[::-1]

    csx = np.cumsum(x)
    csx2 = np.cumsum((x**2).astype("float64"))
    n = np.arange(1, len(x)+1)
    cmeanx = csx/n
    cmx2 = csx2/n
    cstdx = np.sqrt(cmx2-cmeanx**2+1e-3)
    cerrx = cstdx / np.sqrt(n)
    return cmeanx[::-1], cerrx[::-1]