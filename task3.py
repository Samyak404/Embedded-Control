import numpy as np
import matplotlib.pyplot as plt

# Given parameters
D = 0.1  # Pipe Diameter
dz_dp = -100  # Constant pressure gradient
N = 100  # Number of grid points

ap = 0.08  # Temperature expansion coefficient
bt = 5  # Reference temperature coefficient
t0 = 273  # Reference temperature
mu0 = 0.05  # Base viscosity

# Generate radial values
r_values = np.linspace(0, D/2, N)
delta_r = r_values[1] - r_values[0]

# Initialize matrices for the linear system of equations
A = np.zeros((N, N))
b = np.zeros(N)

# Boundary conditions at the center and outer radius of the pipe
A[0, 0] = 1
A[0, 1] = -1
b[0] = 0

A[N-1, N-1] = 1
b[N-1] = 0

# Populate the matrices for the finite difference method
for i in range(1, N-1):
    t = t0 + (bt * (r_values[i]**2))
    mu = mu0 * (1 + (ap * (t - t0)))

    A[i, i] = (-(2 / delta_r**2) + 1 / (r_values[i] * delta_r))
    A[i, i-1] = 1 / delta_r**2 - 1 / (r_values[i] * delta_r)
    A[i, i+1] = (1 / delta_r**2)
    b[i] = dz_dp / mu

# Solve the system of linear equations
u_z_values = np.linalg.solve(A, b)

# Plot the axial velocity profile in the pipe
plt.plot(r_values, u_z_values)
plt.xlabel('Radius (r)')
plt.ylabel('Axial Velocity ($u_z$)')
plt.title('Axial Velocity Profile in a Pipe')
plt.show()

"""
Comments:
1. **Parameter Initialization:** Define parameters such as pipe diameter (`D`), pressure gradient (`dz_dp`), number of grid points (`N`), and coefficients for temperature and viscosity.

2. **Grid Generation:** Create an array of radial values (`r_values`) representing the cross-section of the pipe.

3. **Matrix Initialization:** Set up matrices `A` and `b` to represent the system of linear equations for the finite difference method.

4. **Boundary Conditions:** Implement boundary conditions at the center and outer radius of the pipe.

5. **Finite Difference Method:** Populate the matrices based on the finite difference method, taking into account temperature-dependent viscosity.

6. **Solution:** Solve the system of linear equations using NumPy's `linalg.solve` function.

7. **Plotting:** Plot the axial velocity profile as a function of radial distance.

This code models the flow in a pipe with temperature-dependent viscosity, providing insights into how the axial velocity varies across the pipe's cross-section.
"""