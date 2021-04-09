from ising_model import SimulateMH
import numpy as np
import logging

logger = logging.getLogger()

def to_run(i, steps, T, N, M, freq, SEED,bc=1, return_engine=False, init="random"):
    engine = SimulateMH(N, M, freq, bc, SEED)
    engine.set_T(T)
    if init=="random":
        engine.random_init()
    elif init=="constant":
        engine.constant_init()
    else:
        assert True
    engine.make_steps(steps)
    if return_engine is True:
        return (T, N, engine.get_sampled_M().copy(), engine.get_sampled_E().copy(), engine)
    return (T, N, engine.get_sampled_M().copy(), engine.get_sampled_E().copy())


def findpos(x, start=0):
    def findpos_rec(x, start=0):
        x1 = x[start:]

        mean = np.mean(x1)
        std = np.mean(x1)
        if mean < x[0]:
            mean = - mean
            x = -x
            x1 = -x1
        x2 = x1[x1>=mean]
        mean2 = np.mean(x2)

        #         print(np.where(x>mean2))
        return np.where(x>=mean2)[0][0]
    pos = [0]
    for _ in range(400):
        pos.append(findpos_rec(x, start=pos[-1]))
        if pos[-1]==pos[-2]:
            break
    return pos[-1]

def find_relaxation(T, N,M, steps, SEED):
    try:
        freq=max(1,steps//10**6)
        _,_,Ms, Es, engine = to_run(1, steps, T=T, N=N,M=M, freq=freq,
                                    SEED=SEED, return_engine=True,
                                    init="constant")
        pos = findpos(Ms)
        RTM_const = pos * freq
        pos = findpos(Es)
        RTE_const = pos * freq

        _,_,Ms, Es, engine = to_run(1, steps, T=T, N=N, freq=freq,
                                    SEED=SEED+5, return_engine=True,
                                    init="random")
        pos = findpos(Ms)
        RTM_rand = pos * freq
        pos = findpos(Es)
        RTE_rand = pos * freq
        # return None
        return (T, RTM_const, RTE_const, RTM_rand, RTE_rand)
    except Exception as e:
        logger.error(e)
        return None
