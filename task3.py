import numpy as np
import matplotlib.pyplot as plt

D = 0.1  
dz_dp = -100
N = 100

ap = 0.08
bt = 5
t0 = 273
mu0 = 0.05  

r_values = np.linspace(0, D/2, N)
delta_r = r_values[1] - r_values[0]

t = np.zeros(N)
mu = np.zeros(N)

A = np.zeros((N, N))
b = np.zeros(N)

A[0, 0] = 1
A[0, 1] = -1
b[0] = 0

A[N-1, N-1] = 1
b[N-1] = 0

for i in range(1, N-1):
    t[i] = t0 + (bt * (r_values[i]**2))
    mu[i] = mu0 * (1 + (ap * (t[i] - t0)))

    A[i, i] = (-(2 / delta_r**2) + 1 / (r_values[i] * delta_r))
    A[i, i-1] = 1 / delta_r**2 - 1 / (r_values[i] * delta_r)
    A[i, i+1] = (1 / delta_r**2)
    b[i] = dz_dp / mu[i]

u_z_values = np.linalg.solve(A, b)

plt.plot(r_values, u_z_values)
plt.xlabel('Radius (r)')
plt.ylabel('Axial Velocity ($u_z$)')
plt.title('Axial Velocity Profile in a Pipe')
plt.show()