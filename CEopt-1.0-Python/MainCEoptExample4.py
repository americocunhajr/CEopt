import time
import numpy as np
from scipy.integrate import solve_ivp
from numpy.linalg import norm
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

print(' ------------------- ')
print(' MainCEoptExample4.py ')
print(' ------------------- ')

# Define the Misfit Function
# -----------------------------------------------------------------
def MyMisfitFunc(x, ydata, tspan):
    x = np.atleast_2d(x)    # Ensure x is 2D
    Ns, Nvars = x.shape
    J = np.zeros(Ns)
    wn_arr  = x[:, 0]  # model parameter 1
    ksi_arr = x[:, 1]  # model parameter 2
    y0_arr  = x[:, 2]  # model parameter 3
    v0_arr  = x[:, 3]  # model parameter 4
    
    # Loop over each sample
    for n in range(Ns):
        # Define the ODE system for the nth sample.
        def dydt(t, y, wn, ksi):
            # y is [position, velocity]
            return [y[1], -wn**2 * y[0] - 2 * ksi * wn * y[1]]
        
        # Solve the ODE using solve_ivp.
        sol = solve_ivp(dydt, (tspan[0], tspan[-1]), [y0_arr[n], v0_arr[n]],
                        t_eval=tspan, args=(wn_arr[n], ksi_arr[n]))
        ymodel = sol.y.T  # Transpose so that each row corresponds to a time point.
        J[n] = norm(ydata - ymodel[:, 0]) / np.sqrt(len(ydata))
    return J
# -----------------------------------------------------------------

# System response observations
wn       = 1.0
ksi      = 0.1
y0_true  = 2.0
v0_true  = 0.5
t0       = 0.0
t1       = 30.0
Ndata    = 50
time_vec = np.linspace(t0, t1, Ndata)
w        = wn * np.sqrt(1 - ksi**2)
A        = np.sqrt(y0_true**2 + ((v0_true + ksi * wn * y0_true) / w)**2)
phi      = np.arctan(y0_true * w / (v0_true + ksi * wn * y0_true))
ytrue    = A * np.exp(-ksi * wn * time_vec) * np.sin(w * time_vec + phi)
ydata    = ytrue + 0.05 * y0_true * np.random.randn(Ndata)

# Objective function (minimize the discrepancy)
F = lambda x: MyMisfitFunc(x, ydata, time_vec)

# Bound for design variables
lb = np.array([0.0, 0.0, -3, -3])
ub = np.array([2.0, 1.0,  3,  3])

# Define parameters for the CE optimizer
CEstr = {}
CEstr['Nsamp']  = 50       # Number of samples
CEstr['TolAbs'] = 1.0e-3   # Absolute tolerance
CEstr['TolRel'] = 1.0e-2   # Relative tolerance

# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEobj = CEopt(F, None, None, lb, ub, None, CEstr)
elapsedTime = time.time() - start_time

print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsedTime))
