import numpy as np


def mean_with_err(arr, axis=-1):
    return arr.mean(axis=axis), arr.std(axis=axis) / np.sqrt(arr.shape[axis])


def median_with_err(arr, axis=-1):
    mu = np.median(arr, axis=axis)
    return mu, np.mean(np.abs(arr - mu), axis=axis)
