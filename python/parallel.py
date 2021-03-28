from ising_model import Simulate_MH
import numpy as np


def to_run(i, steps, T, N, freq, SEED, return_engine=False):
    print(i)
    engine = Simulate_MH(int(N), int(N), int(freq), False, int(SEED))
    engine.set_T(T)
    engine.random_init()
    engine.make_steps(int(steps))
    if return_engine is True:
        return (engine.get_sampled_M().copy(), engine.get_sampled_E().copy(), engine)
    return (engine.get_sampled_M().copy(), engine.get_sampled_E().copy())
