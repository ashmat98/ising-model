import numpy as np
from scipy.signal import convolve
from IPython.display import Audio

TTc = 2/np.log(1+np.sqrt(2))

def beep():
    """
    beeps at run. Useful tool in jupyter notebook.
    """
    sound_file = 'http://www.soundjay.com/button/beep-07.wav'
    return Audio(sound_file, autoplay=True)    

def mean_with_err(arr, axis=-1):
    """
    returns mean and error estimate of the array
    """
    return arr.mean(axis=axis), arr.std(axis=axis) / np.sqrt(arr.shape[axis])


def median_with_err(arr, axis=-1):
    """
    returns median and error estimate of the array
    """
    mu = np.median(arr, axis=axis)
    return mu, np.mean(np.abs(arr - mu), axis=axis)


def moving_average(x, w, stride=1):
    """
    Fast moving average with stride.
    """
    s = (x.shape[0] // stride) * stride
    x = x[-s:]
    x = x.reshape((x.shape[0] // stride, stride))
    assert w % stride == 0

    return convolve(x, np.ones((w // stride, stride)), 'valid', method="direct")[:, 0] / w


def moving_mean_err(x, w, stride):
    mean_x_square = moving_average((x ** 2).astype("float64"), w, stride)
    mean_x = moving_average(x, w, stride)
    err = np.sqrt(mean_x_square - mean_x ** 2) / np.sqrt(w)
    return mean_x, err


def cum_mean_err(x):
    x = x[::-1]

    csx = np.cumsum(x)
    csx2 = np.cumsum((x ** 2).astype("float64"))
    n = np.arange(1, len(x) + 1)
    cmeanx = csx / n
    cmx2 = csx2 / n
    cstdx = np.sqrt(cmx2 - cmeanx ** 2 + 1e-3)
    cerrx = cstdx / np.sqrt(n)
    return cmeanx[::-1], cerrx[::-1]


def steps_needed_normalized(T):
    """
    Returns needed steps for relaxation at given temperature and
    for lattice size 1x1. Further multiply with (Nr x Nc)
    """
    steps = 10**5
    b=0.08
    n=0.1
    if T<2.5:
        steps = 3*10**6
    else:
        steps = np.exp(np.log(3*10**6)*(b/(b+T-2.5))**n)
    return steps/32/32
def relaxation_time_normalized(T):
    """
    Returns relaxation time at given temperature and
    for lattice size 1x1. Further multiply with (Nr x Nc)
    """
    a,b,c = 2.99, 1.85, 7.39
    def f(x,a,b,c):
        return a/(x-b)+c
    M = 4.5e5
    Tm = a/(np.log(M)-c) + b
    
    
    if T<Tm:
        steps = M
    else:
        steps = np.exp(f(T, a,b,c)) 
    return steps/32/32


def findpos(x, start=0):
    """
    finds point after which assumed stationary state.
    """
    def findpos_rec(x, start=0):
        x1 = x[start:]

        mean = np.mean(x1)
#         std = np.mean(x1)
        if mean < x[0]:
            mean = - mean
            x = -x
            x1 = -x1
#         mean2 = np.mean(x1[x1>=mean])
        mean2 = mean
        #         print(np.where(x>mean2))
        return np.where(x>=mean2)[0][0]
    pos = [0]
    for _ in range(400):
        pos.append(findpos_rec(x, start=pos[-1]))
        if pos[-1]==pos[-2]:
            break
    return pos[-1]

def calc_autocor(dt, vals):
    if len(vals[dt:]) == 0:
        return 0
    a = np.mean(vals[dt:] * vals[:-dt])
    return a

def find_decorrelation_time(vals):
    """
    Find decorrelation time with binary search.
    """
    mean = np.mean(vals)
    vals1 = vals-mean
    var = np.var(vals)
    
    def check(a):
        return a*np.e > var
    
    last_OK = 0
    add = 1
    while True:
        dt = add + last_OK
        a = calc_autocor(dt, vals1)
        if check(a):
            add *=2
        elif add == 1:
            break
        else:
            last_OK += add//2
            add = 1
    return last_OK

def arrayify(*args):
    """
    converts arguments to numpy arrays
    """
    return (np.array(x) for x in args)