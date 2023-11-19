import numpy as np
import matplotlib.pyplot as plt

# Given parameters
D = 0.1      # Diameter of pipe m
R = D / 2    # Radius of pipe m
mu = 0.01    # Viscosity kg / m.s
dp_dz = -100 # Constant pressure gradient Pa/m
N = 100      # Number of grid points
gamma = 0.25
L = 1.0      # Lenght in z direction

# Discretization
N_r = N
N_z = N

# Generate arrays for radial and axial coordinates
r_values = np.linspace(0, R, N)
z_values = np.linspace(0, L, N)

# Calculate step sizes
delta_r = r_values[1] - r_values[0]
delta_z = z_values[1] - z_values[0]

# Initialize coefficient matrix A and right-hand side vector b
A = np.zeros((N_r * N_z, N_r * N_z))
b = np.zeros(N_r * N_z)

# Loop to populate the coefficient matrix A and vector b
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

# Solve the system of linear equations
u_z = np.linalg.solve(A, b)

# Reshape the solution to a matrix for plotting
u_z_matrix = u_z.reshape((N_r, N_z))

# Create meshgrid for 3D plotting
r, z = np.meshgrid(r_values, z_values)

# Plotting the 3D surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(r, z, u_z_matrix.T, cmap='viridis')
ax.set_xlabel('Radius (m)')
ax.set_ylabel('Axial Position (m)')
ax.set_zlabel('Velocity (m/s)')
ax.set_title('Velocity Profile in a Pipe')
plt.show()

"""
Explanation:

1. **Discretization and Grid Generation:**
   - The radial and axial coordinates are discretized using linspace to create arrays `r_values` and `z_values`.

2. **Coefficient Matrix and Right-hand Side Vector:**
   - A loop is used to populate the coefficient matrix `A` and the right-hand side vector `b` based on the discretized Navier-Stokes equations.

3. **Boundary Conditions:**
   - Special treatment is given to the boundary nodes by setting the diagonal elements of `A` to 1 and corresponding elements in `b` to 0.

4. **Solving the System:**
   - The system of linear equations `Ax = b` is solved using NumPy's `linalg.solve` function.

5. **Reshaping for Plotting:**
   - The solution vector `u_z` is reshaped into a matrix for 3D plotting.

6. **3D Plotting:**
   - Matplotlib is used to create a 3D surface plot of the velocity profile in the pipe.

The provided code generates a 3D plot showing the velocity profile in a pipe considering both radial and axial coordinates. 
"""