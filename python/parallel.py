from ising_model import Simulate_MH
import numpy as np


def to_run(i, steps, T, N, freq, SEED):
    print(i)
    engine = Simulate_MH(int(N), int(N), int(freq), int(SEED))
    engine.set_T(T)
    engine.random_init()
    engine.make_steps(int(steps))
    return (engine.get_sampled_M(), engine.get_sampled_E())
