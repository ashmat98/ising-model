from ising_model import SimulateMH
BC = SimulateMH.BoundaryCondition
import numpy as np
import logging

from utils import *

logger = logging.getLogger()

def to_run(i, steps, T, N, M, freq, SEED,bc=BC.Periodic, return_engine=False, init="random", H=0, omega=0):
    engine = SimulateMH(N, M, freq,H, omega, bc, SEED)
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

def to_run_for_MvsT(i, steps, T, N, M, freq, SEED,bc=BC.Periodic, return_engine=False, init="random", H=0, omega=0):
    mean_ranges = (10**np.linspace(0, np.log10(steps), 20)).astype(int)
    
    T, N, M, E = to_run(i, steps, T, N, M, freq, SEED,bc, return_engine, init, H, omega)
    
    meanMs, errMs = [], []
    meanEs, errEs = [], []
    for mr in mean_ranges:
        mean, err = mean_with_err(M[-mr:])
        meanMs.append(mean), errMs.append(err)
        
        mean, err = mean_with_err(E[-mr:])
        meanEs.append(mean), errEs.append(err)
        

    meanMs, errMs, meanEs, errEs = arrayify(meanMs, errMs, meanEs, errEs)
    return (T, N, meanMs, errMs, meanEs, errEs)



def find_relaxation(T, N, M, steps, SEED):
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

def find_sigma_em(T, N, steps, freq,  SEED, bc=BC.Periodic):
    
    
#     pos1 = findpos(Es)
    pos1 = 0
    pos2 = int(steps_needed_normalized(T)*N*N)
    pos=max(pos1, pos2)
    _,_,Ms, Es = to_run(1, 3*pos + steps, T=T, N=N,M=N, freq=freq,
                                    SEED=SEED+5, return_engine=False,bc=bc,
                                    init="random")
    Es = Es[3*pos//freq:]
    Es = Es.astype("float64")
    if len(Es) == 0:
        return T, len(Es), pos1,pos2, np.nan, np.nan
    
    return (N, T, len(Es), pos1, pos2, np.mean(Es), np.std(Es), np.mean(Es**3), np.mean(Es**4), 
                                       np.mean(Ms), np.std(Ms), np.mean(Ms**3), np.mean(Ms**4))
                                    
find_sigma_em.column_names = ["N", "temp", "len(Es)", "pos1","pos2", "mean_E", "std_E", "E^3", "E^4", 
                                                                     "mean_M", "std_M", "M^3", "M^4"]    

def do_find_decorrelation_time(T, N, M, steps, freq, SEED):
    relax_steps = int(steps_needed_normalized(T)*N*M)
    
    _,_,Ms, Es = to_run(1, 3*relax_steps + steps, T=T, N=N, M=M, freq=freq, bc=BC.Periodic,
                            SEED=SEED, return_engine=False, 
                            init="random")
    pos2 = findpos(Es)*freq
    pos3 = max(pos2,relax_steps)
    
    Ms = Ms[(2*relax_steps)//freq:]
    dcort = find_decorrelation_time(Ms)
    return dcort * freq