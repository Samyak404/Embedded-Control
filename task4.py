import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

D = 0.1 
R = D / 2
mu = 0.01  
dp_dz = -100
N = 100
gamma = 0.25
L = 1.0

N_r = N 
N_z = N

r_values = np.linspace(0, R, N)
z_values = np.linspace(0, L, N)

delta_r = r_values[1] - r_values[0]
delta_z = z_values[1] - z_values[0]

A = np.zeros((N_r * N_z, N_r * N_z))
b = np.zeros(N_r * N_z)


for i in range(N_r):
    for j in range(N_z):
        idx = i + j * N_r

        if i == 0 or i == N_r - 1 or j == 0 or j == N_z - 1:
            A[idx, idx] = 1
            b[idx] = 0
        else:
            A[idx, idx + N_z] = mu / (delta_r ** 2) + mu / (delta_r * r_values[i]) + gamma / (delta_z ** 2)
            A[idx, idx] = -2 * mu / (delta_r ** 2) - mu / (delta_r * r_values[i]) - 2 * gamma / (delta_z ** 2)
            A[idx, idx - N_z] = mu / (delta_r ** 2) + gamma / (delta_z ** 2)
            b[idx] = dp_dz


u_z = np.linalg.solve(A, b)

u_z_matrix = u_z.reshape((N_r, N_z))

print(u_z_matrix.shape)

r, z = np.meshgrid(r_values, z_values)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(r, z, u_z_matrix.T, cmap='viridis')
ax.set_xlabel('Radius (m)')
ax.set_ylabel('Axial Position (m)')
ax.set_zlabel('Velocity (m/s)')
ax.set_title('Velocity Profile in a Pipe')
plt.show()