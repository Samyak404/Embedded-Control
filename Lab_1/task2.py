# Import the necessary libraries

import numpy as np
# Numpy provides functions to solve matrices

import matplotlib.pyplot as plt
# Mathplotlib provides functions to plots graphs

# Initialize variables

D = 0.1 
# Pipe Diameter in meters

mu = 0.01
# Viscosity of the fluid in kg/ms

dp_dz = -100
# Pressure gradient in Pa/m

N = 100
# Number of grid points

# Generate radial grid points by discretizing

r_values = np.linspace(0, D/2, N)
# Create N descrete partitions from 0 to R (D/2)

delta_r = r_values[1] - r_values[0]
# Difference between 2 adjecent grid points

# Initialize matrix A and vector b
A = np.zeros((N, N))
b = np.zeros(N)

# Boundary conditions

# at r = 0
A[0, 0] = 1
A[0, 1] = -1
b[0] = 0

# at r = R (D/2)
A[N-1, N-1] = 1
b[N-1] = 0

# Using finite difference method fill in matrix A and vector b
for i in range(1, N-1):
    A[i, i] = (-(2 / delta_r**2) + 1 / (r_values[i] * delta_r))
    A[i, i-1] = 1 / delta_r**2 - 1 / (r_values[i] * delta_r)
    A[i, i+1] = (1 / delta_r**2)
    b[i] = dp_dz / mu

# Solving linear equations using numpy function
u_z_values = np.linalg.solve(A, b)

# Plot the velocity profile
plt.plot(r_values, u_z_values)
plt.xlabel('Radius (r)')
plt.ylabel('Axial Velocity ($u_z$)')
plt.title('Axial Velocity Profile in a Pipe')
plt.show()