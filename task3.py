import numpy as np
import matplotlib.pyplot as plt

D = 0.1
N = 10
dz_dp = -100 

t0 = 273
ap = 0.08
bt = 5
mu0 = 0.05

r_values = np.linspace(0, D/2, N) 
delta_r = r_values[2] - r_values[1]

A = np.zeros((N, N))
b = np.zeros(N)

for i in range(1, N-1):
    t = t0 + bt * r_values[i]**2
    mu = mu0 * (1 + ap * (t - t0))

    A[i, i] = (-2*mu) / (delta_r**2) + mu / ((i+1) * delta_r)
    A[i, i-1] = mu / (delta_r**2) + mu / ((i+1) * delta_r)
    A[i, i+1] = mu / (delta_r**2) + mu / ((i+1) * delta_r)
    b[i] = dz_dp

A[0, 0] = 1
A[N-1, N-1] = 1
b[0] = 15 
b[N-1] = 0  

u_z_values = np.linalg.solve(A, b)

plt.plot(r_values, u_z_values)
plt.show()
