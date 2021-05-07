# Ising-model
This engine runs Metropolis-Hasting algorithm on 2D ising model. It provides base framework for the extension to more 
complex models with extended interactions. At this stage, only 4-neighbour interactions are accounted.

The project strcture is:
```
ising-model
 ┣ cpp                      /* C++ part, simulation engine */
 ┃ ┣ base.cpp 
 ┃ ┣ base.h                 /* base class of the engine */
 ┃ ┣ CMakeLists.txt         /* configuration for CMake */
 ┃ ┣ ising-model.cpp        /* Language bindings using PyBind11 */
 ┃ ┣ main.cpp               /* sample usage of the engine inside C++ */
 ┃ ┣ metropolis-hasting.cpp 
 ┃ ┗ metropolis-hasting.h   /* Metropolis-Hasting algorithm */
 ┣ python
 ┃ ┣ dumps                  /* simulated data, which is not uploaded to 
 ┃ ┃                            the github */
 ┃ ┃ ┗ ...
 ┃ ┣ figs                   /* generated figures */
 ┃ ┃ ┗ ...
 ┃ ┣ parallel.py            /* functions that run the engine for specific
 ┃ ┃                            task, can be run in paralel. */
 ┃ ┣ style.py               /* style for plots */
 ┃ ┗ utils.py               /* utility functions */
 ┣ CMakeLists.txt           /* configuration for CMake */
 ┗ README.md                /* This file */
```

###Compilation
1. Install [conda environment](https://docs.conda.io/en/latest/miniconda.html) for convenience.
2. Run these commands in the shell
```shell
$ conda install pybind11 cmake
$ git clone ...
$ cd ising-model
$ mkdir build
$ cd build
$ cmake cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```
3. This will create the python module. E.g. on Windows something like 
   `ising_model.cp38-win_amd64.pyd` or e.g. `ising_model.cpython-37m-x86_64-linux-gnu.so` on Linux.
   This can be imported to Python (see below).
### Usage
Sample usage in `C++`
```c++
SimulateMH engine(32, 32, 1, 1);
engine.set_T(10);
engine.make_steps(stps);
auto E = engine.get_sampled_E();
auto M = engine.get_sampled_M();
etc.
```

Sample usage in `Python`:
```python
from ising_model import SimulateMH
BC = SimulateMH.BoundaryCondition
engine = SimulateMH(Nr=32, Mc=32, frequency_to_store=1, H=1, omega=0, 
                    bc=BC.Periodic, SEED=42)
engine.random_init()
# or engine.constant_init()
engine.make_steps(10**8) # ~10 seconds runtime
M = engine.get_sampled_M()
E = engine.get_sampled_E()
etc.
```

Paralelisation usnig `multiprocessing.pool`:
```python
from parallel import to_run

```

