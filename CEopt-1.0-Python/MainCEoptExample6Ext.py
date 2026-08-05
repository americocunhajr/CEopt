# -----------------------------------------------------------------
#  MainCEoptExample6Ext.py
# -----------------------------------------------------------------
#  programmer: Americo Cunha
#              americo.cunhajr@gmail.com

# 
#  Originally programmed in: Mar 28, 2025
#           Last updated in: jun 12, 2025
# -----------------------------------------------------------------
#  Example 6: Nonsmooth function with conic constraints
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
#np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample6Ext.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
# Objetive function
# -----------------------------------------------------------------
def PatternSearchFunc(x):

    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]

    nSamples = x.shape[0]
    F        = np.zeros((nSamples, 1))

    for i in range(nSamples):
        if x1[i] < -5:
            F[i,0] = (x1[i] + 5)**2 + np.abs(x2[i])
        elif x1[i] < -3:
            F[i,0] = -2 * np.sin(x1[i]) + np.abs(x2[i])
        elif x1[i] < 0:
            F[i,0] = 0.5 * x1[i] + 2 + np.abs(x2[i])
        else:  # x1[i] >= 0
            F[i,0] = 0.3 * np.sqrt(x1[i]) + 5/2 + np.abs(x2[i])

    return F
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# Constraint function
# -----------------------------------------------------------------
def ConicConstraints(x):

    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]

    nSamples = x.shape[0]
    G        = np.zeros((nSamples, 1))
    H        = np.zeros((nSamples, 1))

    for i in range(nSamples):
        G[i,0] = 2.0 * x1[i]**2 + x2[i]**2 - 3.0
        H[i,0] = (x1[i] + 1.0)**2 - (x2[i] / 2.0)**4
    
    return G, H
# -----------------------------------------------------------------

# Bound for design variables and initial mean
lb  = np.array([-6, -4])
ub  = np.array([ 2,  4])
mu0 = np.array([-4,  2])

# Objective function and constraints
F       = PatternSearchFunc
nonlcon = ConicConstraints

# Cross-entropy optimizer struct
CEstr = {}
CEstr['isVectorized'] = True    # vectorized function
CEstr['TolCon'      ] = 1.0e-6  # constraint tolerance

# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(F, mu0, None, lb, ub, nonlcon, CEstr)
elapsedTime = time.time() - start_time

print("\n")
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsedTime))

# Check constraint violation
G, H = nonlcon(Xopt)
print("\n")
print("--------------------------------------------------------")
print("Check if the inequality constraints is <= 0")
print("--------------------------------------------------------")
print("G(x) = {:+.6E}".format(G[0,0]))
print("--------------------------------------------------------")

print("\n")
print("--------------------------------------------------------")
print("Check if the equality constraints is = 0")
print("--------------------------------------------------------")
print("H(x) = {:+.6E}".format(H[0,0]))
print("--------------------------------------------------------")

# Meshgrid for visualization
Ngrid  = 100
xRange = np.linspace(lb[0], ub[0], Ngrid)
yRange = np.linspace(lb[1], ub[1], Ngrid)
X, Y   = np.meshgrid(xRange, yRange)
pts    = np.column_stack((X.ravel(), Y.ravel()))
Z      = F(pts)
Z      = Z.reshape(X.shape)
C, Ceq = nonlcon(pts)
C      = C.reshape(X.shape)
Ceq    = Ceq.reshape(X.shape)

# Custom colors and colormaps
MyYellow    = [0.9290, 0.6940, 0.1250]
MyBlue      = [0.0000, 0.4470, 0.7410]
MyRed       = [0.6350, 0.0780, 0.1840]
Ncolors     = 50
MyColorMap1 = colormaps['viridis'].resampled(Ncolors)
MyColorMap2 = colormaps['inferno'].resampled(Ncolors)

# Set up the figure for visualization
plt.figure()
plt.clf()

# Plot the objective function as filled contours.
cont1 = plt.contourf(X, Y, Z, levels=Ncolors, cmap=MyColorMap1)
plt.colorbar()

# Plot the equality constraint H(x)=0 as a contour line.
cont2 = plt.contour(X, Y, Ceq, levels=[0], colors='k', linewidths=2)

# Plot the inequality constraint G(x) ≤ 0.
# Create a mask for feasible region: G(x) ≤ 0
mask = (C <= 0)

# Create a filled array (1 inside feasible region, NaN outside)
filled = np.where(mask, 1, np.nan)

# Use contourf to fill the region between 0.5 and 1.5
cont3 = plt.contourf(X, Y, filled, levels=[0.5, 1.5], colors=[MyYellow], alpha=0.7)

# Plot the initial approximation.
fig4, = plt.plot(mu0[0], mu0[1], 'D', color=MyRed, markersize=7,
                 markerfacecolor=MyRed, label='Initial Approximation')
# Plot CE iterations.
for j in range(CEstr['iter']):
    plt.plot(CEstr['xmean'][j, 0], CEstr['xmean'][j, 1], 'D', color=MyRed,
             markersize=7, label='CE Iterations' if j == 0 else "_nolegend_")

# Plot the CE solution.
fig6, = plt.plot(Xopt[0, 0], Xopt[0, 1], 'x', color=MyBlue, linewidth=3,
                 markersize=25, label='CE Optimum')

# Labeling
plt.xlabel('x1', fontsize=20, fontname='Arial')
plt.ylabel('x2', fontsize=20, fontname='Arial')

# Create proxy artists for legend
from matplotlib.lines import Line2D
proxy_eq = Line2D([0], [0], color='k', lw=2, label='Equality Constraint')
proxy_in = Line2D([0], [0], marker='s', color=MyYellow, markersize=15,
                  label='Inequality Constraint', linestyle='None')
proxy_init = Line2D([0], [0], marker='D', color=MyRed, markersize=7,
                    label='Initial Approximation', linestyle='None')
proxy_iter = Line2D([0], [0], marker='D', color=MyRed, markersize=7,
                    label='CE Iterations', linestyle='None')
proxy_sol = Line2D([0], [0], marker='x', color=MyBlue, markersize=25,
                   label='CE Optimum', linestyle='None')
plt.legend(handles=[proxy_eq, proxy_in, proxy_init, proxy_iter, proxy_sol],
           loc='lower left', fontsize=12)

# Set font for tick labels
ax = plt.gca()
for label in ax.get_xticklabels():
    label.set_fontname('Arial')
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
ax.tick_params(labelsize=18)

plt.box(True)
plt.tight_layout()

# Save the figure as PNG
plt.savefig('CEoptExample6.png', format='png', dpi=300)

# Show the figure
plt.show()
