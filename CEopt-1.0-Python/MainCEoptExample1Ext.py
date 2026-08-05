# -----------------------------------------------------------------
#  MainCEoptExample1Ext.py
# -----------------------------------------------------------------
#  programmer: Americo Cunha
#              americo.cunhajr@gmail.com
# 
#  Originally programmed in: Mar 27, 2025
#           Last updated in: Mar 28, 2025
# -----------------------------------------------------------------
#  Example 1: A Gaussian mixture in 1D
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from CEopt import CEopt  # Assumes CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample1Ext.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
# Objective function
# -----------------------------------------------------------------
def F(x):
    # x can be a scalar or a NumPy array; operations are elementwise.
    return -0.8 * np.exp(-(x - 2)**2) - 0.5 * np.exp(-(x + 2)**2) + 1
# -----------------------------------------------------------------

# Bound for design variables
lb = -5
ub =  5

# Initialize mean and std. dev. vectors
mu0    = (ub + lb) / 2.0
sigma0 = (ub - lb) / 6.0

# Run the CE optimizer and time the execution
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(F, mu0, sigma0, lb, ub)
elapsed_time = time.time() - start_time

print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsed_time))

# Domain for plotting
xgrid = np.linspace(lb, ub, 1000)
ygrid = F(xgrid)

# Color map: use 15 discrete colors from the 'magma' colormap
Ncolors = 15
cmap    = colormaps['magma'].resampled(Ncolors)

# Set up the figure
plt.figure()
plt.clf()
plt.rcParams['font.size'] = 18

# Plotting the Gaussian distributions from the CE history.
# In Matlab the loop is "for i = 1:4:N" (1-indexed); in Python we use 0-indexing.
for i in range(0, Ncolors, 4):
    # For a one-variable problem, CEstr['xmean'] and CEstr['sigma'] have shape (iterations, 1)
    mu    = CEstr['xmean'][i, 0]
    sigma = CEstr['sigma'][i, 0]
    # Compute the normalized Gaussian PDF with mean mu and std sigma
    y_gauss = np.exp(- (xgrid - mu)**2 / (2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)
    # Select color: mimic Matlab’s cmap(N+1-i,:) (with 1-indexing)
    color = cmap(Ncolors - 1 - i)
    plt.plot(xgrid, y_gauss, color=color, linewidth=0.8)

# Plot the objective function F(x)
plt.plot(xgrid, ygrid, 'b-', linewidth=2)

# Highlight the optimal point with a black filled circle
plt.plot(Xopt, Fopt, 'ko', markerfacecolor='k', markersize=8)

# Labeling
plt.xlabel('x', fontsize=20)
plt.ylabel('F(x)', fontsize=20)

# Setting plot limits
plt.xlim([lb, ub])
plt.ylim([0, 25])

plt.box(True)
plt.tight_layout()

# Save the figure as PNG
plt.savefig('CEoptExample1.png', format='png', dpi=300)

# Show the figure
plt.show()
