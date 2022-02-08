# Import packages 
using Plots
using LinearAlgebra
gr()

# User variables 
Nx = 50;    # Number of spatial grid points
L = 1.0;    # Length of 1D space

u0 = 10.0;  # LHS boundary condition
uL = 0;     # RHS boundary condition

Δt = 0.1;   # Timestep size
t = 20.0;   # Time to advance to

K = 2.0;    # Heat coeffcient

# Set up 
x = range(0,length=Nx,stop=L)   # Generate x-axis gridding
Nt = t/Δt;     # Calculate number of time steps needed
u_init = zeros(Nx-2,1); # Define initial temp values

# Set up operator matrix (blunt method)
L = zeros(Nx-2,Nx-2);
for i in 1:(Nx-2)
    if (i == 1)
        L[1,1] = -2;
        L[1,2] = 1;
    elseif (i == (Nx-2))
        L[(Nx-2),(Nx-2)] = -2;
        L[(Nx-2),(Nx-2)-1] = 1;
    else
        L[i,i] = -2;
        L[i,(i+1)] = 1;
        L[i,(i-1)] = 1;
    end
end

# Set up BC matrix
g = zeros(Nx-2,1);  # Create a matrix of zeros
g[1,1] = u0;        # Add the LHS bc 
g[Nx-2,1] = uL;     # Add the RHS bc 

# Timestepping
# Set up matrices
f = 0.5*K*Δt;  # Calculate factor in front (Kt/2)
A = (I - f*L)       # Calculate A
b = (I + f*L)*u_init + f*2*g  # Calculate b intially

# Solve and advance
u_temp = A\b;  # First timestep, solve for u
for t in 2:Nt  # Run through all timesteps
    global b = (I + f*L)*u_temp + f*2*g  # Update b matrix
    global u_temp = A\b;  # Update u values
end

# Construct final u vector
u_final = zeros(Nx,1)
u_final[2:Nx-1,1] = u_temp
u_final[1,1] = u0;
u_final[Nx,1] = uL;

# Plot final result
plot(x, u_final)
xaxis!("x position")
yaxis!("Temp")


