from ising_model import Simulate_MH
import numpy as np


def to_run(steps):
    engine = Simulate_MH(32, 32, 100, np.random.randint(0, 1000000))
    engine.set_T(100)
    engine.random_init()
    engine.make_steps(steps)
    return (engine.get_sampled_M(), engine.get_sampled_E())
