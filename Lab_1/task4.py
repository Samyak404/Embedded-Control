import numpy as np
import matplotlib.pyplot as plt

# Initialize variables
N = 100    # Number of grid points
N_r = N    # Along radial direction
N_z = N    # Along axial direction
D = 0.1    # Pipe Diameter in meters
R = D/2    # Radius of pipe Radius = Diameter/2
L = 1.0    # Lenght of pipe in meter
dpdz = -100    # Pressure gradient in Pa/m

# Generate arrays for radial and axial coordinates
r_values = np.linspace(0, R, N_r)
z_values = np.linspace(0, L, N_z)

# Calculate step sizes
delta_r = r_values[1] - r_values[0]
delta_z = z_values[1] - z_values[0]

# Initialize matrix A and vector b
A = np.zeros((N_r * N_z, N_r * N_z))
b = np.full((N_r * N_z),dpdz)

reversed_u_z = np.zeros([N_r,N_z])

mu = 0.01   # Viscosity of the fluid in kg/ms
gamma = 0.25    # A constant

def Index_Function(a,b):
    index = a + b *N_r
    return index

def Reverse_Index(index):
    Y_Index = int(index / N_r)
    X_Index = int(index % N_r)
    return X_Index,Y_Index
    
# Populate the elements of matrix A and vector b
for j in range(N_z):
    for i in range(N_r):
        radius = i * delta_r
        if j == 0:  # Boundary condition at z = 0
            A[Index_Function(i,j), Index_Function(j,i)] = 1
            b[Index_Function(i,j)] = 15

        elif j == (N_z - 1):  # Boundary condition at z = L
            A[Index_Function(i,j), Index_Function(j,i)] = 1
            b[Index_Function(i,j)] = 15

        elif i == 0:  # Boundary condition at r = 0 (axis)
            A[Index_Function(i,j), Index_Function(j,i+1)] = 1
            A[Index_Function(i,j), Index_Function(j,i)] = -1
            b[Index_Function(i,j)] = 0

        elif i == (N_r - 1):  # Boundary condition at r = R (disk)
            A[Index_Function(i,j), Index_Function(j,i)] = 1
            b[Index_Function(i,j)] = 0

        else:
            A[Index_Function(i,j), Index_Function(j,i+1)] = ((mu) / (delta_r ** 2 )) + ((mu) / (radius * delta_r))
            A[Index_Function(i,j), Index_Function(j,i)] = -1 * (((2 * mu) / (delta_r ** 2)) + (mu / (radius * delta_r)) + (2 * gamma / (delta_z ** 2)))
            A[Index_Function(i,j), Index_Function(j,i-1)] = mu / (delta_r ** 2)
            A[Index_Function(i,j), Index_Function(j+1,i)] = gamma / (delta_z ** 2)
            A[Index_Function(i,j), Index_Function(j-1,i)] = gamma / (delta_z ** 2)

# Solving linear equations using numpy function
u_z = np.linalg.solve(A, b)

# Reshape velocity array for plotting
for a, b in enumerate (u_z):
    reversed_u_z_Index = Reverse_Index(a)
    reversed_u_z[reversed_u_z_Index[0]][reversed_u_z_Index[1]] = b

# Create meshgrid for 3D plotting
r, z = np.meshgrid(r_values, z_values)

# Plot the 3D surface plot
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(z, r, reversed_u_z, cmap='viridis')

ax.set_title('3D Velocity Profile in a Pipe')
ax.set_xlabel('Axial Position (m)')
ax.set_ylabel('Radius (m)')
ax.set_zlabel('Velocity (m/s)')

# Set the view angle
ax.view_init(elev=30, azim=45)

# color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()