import numpy as np
import matplotlib.pyplot as plt

# Define parameters
D = 0.1    # Pipe Diameter in meters
mu = 0.01  # Viscosity of the fluid in kg / (m * s)
dz_dp = -100  # Constant pressure gradient in Pa / m
N = 100    # Number of grid points

# Generate radial grid points
r_values = np.linspace(0, D/2, N)
delta_r = r_values[1] - r_values[0]

# Initialize coefficient matrix A and right-hand side vector b
A = np.zeros((N, N))
b = np.zeros(N)

# Boundary conditions at r=0
A[0, 0] = 1
A[0, 1] = -1
b[0] = 0

# Boundary conditions at r=D/2
A[N-1, N-1] = 1
b[N-1] = 0

# Fill the coefficient matrix A and right-hand side vector b using finite difference method
for i in range(1, N-1):
    A[i, i] = (-(2 / delta_r**2) + 1 / (r_values[i] * delta_r))
    A[i, i-1] = 1 / delta_r**2 - 1 / (r_values[i] * delta_r)
    A[i, i+1] = (1 / delta_r**2)
    b[i] = dz_dp / mu

# Solve the system of linear equations
u_z_values = np.linalg.solve(A, b)

# Plot the axial velocity profile
plt.plot(r_values, u_z_values)
plt.xlabel('Radius (r)')
plt.ylabel('Axial Velocity ($u_z$)')
plt.title('Axial Velocity Profile in a Pipe')
plt.show()

''' Explanation of Code
1. **Parameters**: The necessary parameters such as pipe diameter (`D`), fluid viscosity (`mu`), constant pressure gradient (`dz_dp`), and the number of grid points (`N`).
2. **Grid Generation**: Generate radial grid points (`r_values`) using `np.linspace` based on the pipe diameter and the number of grid points.
3. **Coefficient Matrix and Vector Initialization**: Initialize the coefficient matrix `A` and the right-hand side vector `b` with zeros.
4. **Boundary Conditions**: Set up the boundary conditions at `r=0` and `r=D/2` in the coefficient matrix `A` and vector `b`.
5. **Finite Difference Method**: Use a loop to fill the coefficient matrix `A` and vector `b` using finite difference method based on the discretized Navier-Stokes equation.
6. **System Solution**: Solve the system of linear equations `A u = b` using `np.linalg.solve`.
7. **Plotting**: Plot the axial velocity profile as a function of radial distance using `matplotlib`.
This code models the steady-state and fully-developed axial flow of an incompressible Newtonian fluid in a cylindrical pipe and visualizes the axial velocity profile.
'''