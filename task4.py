# Import the necessary libraries

import numpy as np
# Numpy provides functions to solve matrices

import matplotlib.pyplot as plt
# Mathplotlib provides functions to plots graphs

# Initialize variables

D = 0.1 
# Pipe Diameter in meters

R = D/2
# Radius of pipe Radius = Diameter/2

mu = 0.01
# Viscosity of the fluid in kg/ms

dp_dz = -100
# Pressure gradient in Pa/m

N = 100
# Number of grid points

gamma = 0.25
# A constant

L = 1.0
# Lenght og pipe in meter

# Discretization

N_r = N
# Along radial direction

N_z = N
# Along axial direction

# Generate arrays for radial and axial coordinates

r_values = np.linspace(0, R, N)
# Create N descrete partitions from 0 to R along radius

z_values = np.linspace(0, L, N)
# Create N descrete partitions from 0 to L axially

# Calculate step sizes

delta_r = r_values[1] - r_values[0]
# Difference between 2 adjecent grid points along radius

delta_z = z_values[1] - z_values[0]
# Difference between 2 adjecent grid points axially

# Initialize matrix A and vector b
A = np.zeros((N_r * N_z, N_r * N_z))
b = np.zeros(N_r * N_z)

# Populate the elements of matrix A and vector b
for i in range(N_r):
    for j in range(N_z):
        idx = i + j * N_r

        # Boundary conditions
        if i == 0 or i == N_r - 1 or j == 0 or j == N_z - 1:
            A[idx, idx] = 1
            b[idx] = 0
        else:
            A[idx, idx + N_z] = mu / (delta_r ** 2) + mu / (delta_r * r_values[i]) + gamma / (delta_z ** 2)
            A[idx, idx] = -2 * mu / (delta_r ** 2) - mu / (delta_r * r_values[i]) - 2 * gamma / (delta_z ** 2)
            A[idx, idx - N_z] = mu / (delta_r ** 2) + gamma / (delta_z ** 2)
            b[idx] = dp_dz

# Solving linear equations using numpy function
u_z = np.linalg.solve(A, b)

# Reshape the solution to a matrix for plotting
u_z_matrix = u_z.reshape((N_r, N_z))

# Create meshgrid for 3D plotting
r, z = np.meshgrid(r_values, z_values)

# Plot the 3D surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(r, z, u_z_matrix.T, cmap='viridis')
ax.set_xlabel('Radius (m)')
ax.set_ylabel('Axial Position (m)')
ax.set_zlabel('Velocity (m/s)')
ax.set_title('Velocity Profile in a Pipe')
plt.show()