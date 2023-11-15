import numpy as np
import matplotlib.pyplot as plt

D = 0.1
N = 100
mu = 0.01
dz_dp = -100
L = 1
gamma = 0.25

r_values = np.linspace(0, D/2, N)
delta_r = r_values[1] - r_values[0]
delta_z = L/N

A = np.zeros((N*N, N*N))
b = np.zeros(N*N)

for i in range(N):
    for j in range(N):
        index = i * N + j
        r = r_values[i]
        z = j * delta_z
        if i == 0:  # left boundary
            A[index, index] = 1
            b[index] = 15
        elif i == N-1:  # right boundary
            A[index, index] = 1
            b[index] = 15
        elif j == 0:  # bottom boundary
            A[index, index] = -2/(delta_r**2) - 1/(r*delta_r)
            A[index, index+N] = 1/(delta_z**2)
            if r != 0:  # to avoid division by 0
                A[index, index+1] = 1/(2*r*delta_r)
                A[index, index-1] = -1/(2*r*delta_r)
            b[index] = dz_dp/mu - gamma/(delta_z**2)*z
        elif j == N-1:  # top boundary
            A[index, index] = -2/(delta_r**2) - 1/(r*delta_r)
            A[index, index-N] = 1/(delta_z**2)
            if r != 0:
                A[index, index+1] = 1/(2*r*delta_r)
                A[index, index-1] = -1/(2*r*delta_r)
            b[index] = dz_dp/mu - gamma/(delta_z**2)*z
        else:  # interior nodes
            A[index, index] = -2/(delta_r**2) - 2/(delta_z**2)
            A[index, index+1] = 1/(2*r*delta_r)
            A[index, index-1] = -1/(2*r*delta_r)
            A[index, index+N] = 1/delta_z**2
            A[index, index-N] = 1/delta_z**2

u_z_values = np.linalg.solve(A, b)
u_z_values = np.reshape(u_z_values, (N,N))

r_values  = np.linspace(0, D/2, N)
z_values  = np.linspace(0, L, N)

R, Z = np.meshgrid(r_values, z_values)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(R, Z, u_z_values, cmap='coolwarm')
ax.set_xlabel('Radius (m)')
ax.set_ylabel('Axial Position (m)')
ax.set_zlabel('Axial Velocity ($u_z$)')
ax.set_title('3D plot of Velocity Profile in a Pipe')

plt.show()
