from ising_model import SimulateMH
import numpy as np
import logging

from utils import *

logger = logging.getLogger()

def to_run(i, steps, T, N, M, freq, SEED,bc=1, return_engine=False, init="random"):
    engine = SimulateMH(N, M, freq, bc, SEED)
    engine.set_T(T)
    if init=="random":
        engine.random_init()
    elif init=="constant":
        engine.constant_init()
    else:
        assert False
    engine.make_steps(steps)
    if return_engine is True:
        return (T, N, engine.get_sampled_M().copy(), engine.get_sampled_E().copy(), engine)
    return (T, N, engine.get_sampled_M().copy(), engine.get_sampled_E().copy())




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

        _,_,Ms, Es, engine = to_run(1, steps, T=T, N=N,M=M, freq=freq,
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

def find_sigma_e(T, N, steps, freq, SEED):
    
    _,_,Ms, Es = to_run(1, steps, T=T, N=N,M=N, freq=freq,
                                    SEED=SEED+5, return_engine=False,bc=1,
                                    init="random")
    pos1 = findpos(Es)
    pos2 = int(steps_needed_normalized(T)*N*N)
    pos=max(pos1, pos2)
    Es = Es[3*pos:]
    if len(Es) == 0:
        return T, len(Es), pos1,pos2, np.nan, np.nan
    
    return N, T, len(Es), pos1, pos2, np.mean(Es), np.std(Es)
    

def do_find_decorrelation_time(T, N, M, steps, freq, SEED):
    relax_steps = int(steps_needed_normalized(T)*N*M)
    
    _,_,Ms, Es = to_run(1, 3*relax_steps + steps, T=T, N=N, M=M, freq=freq, bc=1,
                            SEED=SEED, return_engine=False, 
                            init="random")
    pos2 = findpos(Es)*freq
    pos3 = max(pos2,relax_steps)
    
    Ms = Ms[(2*relax_steps)//freq:]
    dcort = find_decorrelation_time(Ms)
    return dcort * freq