import ising_model
import numpy as np

import ising_model

for N in range(89,1,-1):
	print("start ",N)
	engine = ising_model.SimulateMH(N,N,1,1,np.random.randint(0,100000))
	engine.set_T(10)
	engine.make_steps(1000)

exit(0)
import numpy as np




eng = ising_model.Simulate_MH(4, 5)

R = np.random.randint(0, 2, size=(eng.Nr, eng.Nc))
eng.set_state(R)
print(R)
print()

for i in range(0, eng.Nr + 2):
    for j in range(0, eng.Nc + 2):
        print(eng.get(i, j), end="")
    print()
print()

x = eng.get_state()
assert np.allclose(x, R)
print(x)
print(x.dtype)
