# -----------------------------------------------------------------
#  MainCEoptExample4Ext.py
# -----------------------------------------------------------------
#  programmer: Americo Cunha
#              americo.cunhajr@gmail.com
# 
#  Originally programmed in: Mar 27, 2025
#           Last updated in: Mar 28, 2025
# -----------------------------------------------------------------
#  Example 4: Identification of a harmonic oscillator
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numpy.linalg import norm
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample4Ext.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
# Misfit Function
# -----------------------------------------------------------------
def MyMisfitFunc(x, ydata, tspan):
    x = np.atleast_2d(x)
    Ns, Nvars = x.shape
    J = np.zeros(Ns)
    wn_arr  = x[:, 0]
    ksi_arr = x[:, 1]
    y0_arr  = x[:, 2]
    v0_arr  = x[:, 3]
    
    for n in range(Ns):
        def dydt(t, y, wn, ksi):
            return [y[1], -wn**2 * y[0] - 2 * ksi * wn * y[1]]
        sol = solve_ivp(dydt, (tspan[0], tspan[-1]), [y0_arr[n], v0_arr[n]],
                        t_eval=tspan, args=(wn_arr[n], ksi_arr[n]))
        ymodel = sol.y.T  # Each row corresponds to a time point.
        J[n] = norm(ydata - ymodel[:, 0]) / np.sqrt(len(ydata))
    return J
# -----------------------------------------------------------------

# System response observations
wn_true  = 1.0
ksi_true = 0.1
y0_true  = 2.0
v0_true  = 0.5
t0       = 0.0
t1       = 30.0
Ndata    = 50
time_vec = np.linspace(t0, t1, Ndata)
w        = wn_true * np.sqrt(1 - ksi_true**2)
A        = np.sqrt(y0_true**2 + ((v0_true + ksi_true*wn_true*y0_true)/w)**2)
phi      = np.arctan(y0_true*w/(v0_true+ksi_true*wn_true*y0_true))
ytrue    = A * np.exp(-ksi_true*wn_true*time_vec) * np.sin(w*time_vec+phi)
ydata    = ytrue + 0.05*y0_true * np.random.randn(Ndata)

# Objective function (minimize discrepancy)
F = lambda x: MyMisfitFunc(x, ydata, time_vec)

# Bound for design variables (4 parameters)
lb = np.array([0.0, 0.0, -3, -3])
ub = np.array([2.0, 1.0,  3,  3])

# CE optimizer parameters
CEstr = {}
CEstr['Nsamp' ] = 50       # Number of samples
CEstr['TolAbs'] = 1.0e-3   # Absolute tolerance
CEstr['TolRel'] = 1.0e-2   # Relative tolerance

# Run the CE optimizer; passing None for initial mean and sigma so defaults are computed.
start_time = time.time()
Xopt, Fopt, ExitFlag, CEobj = CEopt(F, None, None, lb, ub, None, CEstr)
elapsedTime = time.time() - start_time

print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsedTime))

# -----------------------------------------------------------------
# Animation of the parameter identification process
# -----------------------------------------------------------------
for n in range(CEobj['iter']):  # Matlab loop: for n=1:CEobj.iter (Python: 0-indexed)
    # Extract model parameters from the current iteration
    wn_est  = CEobj['xmean'][n, 0]
    ksi_est = CEobj['xmean'][n, 1]
    y0_est  = CEobj['xmean'][n, 2]
    v0_est  = CEobj['xmean'][n, 3]
    
    # Define the ODE model for the estimated parameters
    def dydt(t, y, wn, ksi):
        return [y[1], -wn**2 * y[0] - 2 * ksi * wn * y[1]]
    
    IC = [y0_est, v0_est]
    sol = solve_ivp(lambda t, y: dydt(t, y, wn_est, ksi_est),
                    (time_vec[0], time_vec[-1]), IC, t_eval=time_vec)
    ymodel = sol.y.T
    y_model = ymodel[:, 0]
    
    plt.figure(1)
    plt.clf()
    plt.plot(time_vec, ytrue, '--k', linewidth=2.0, label='real system')
    plt.plot(time_vec, ydata, 'or', linewidth=2.0, label='noisy data')
    plt.plot(time_vec, y_model, '-b', linewidth=2.0, label='identified model')
    plt.ylim([-3, 3])
    plt.legend(fontsize=16)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    # Save the current figure; here we save as PNG
    saveFileName = 'CEoptExample4_{}.png'.format(n+1)
    plt.savefig(saveFileName, dpi=300)
    plt.close(1)  # Close figure 1 to avoid too many open figures
# -----------------------------------------------------------------
