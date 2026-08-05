# -----------------------------------------------------------------
#  MainCEoptExample3Ext.py
# -----------------------------------------------------------------
#  programmer: Americo Cunha
#              americo.cunhajr@gmail.com
# 
#  Originally programmed in: Mar 27, 2025
#           Last updated in: Mar 28, 2025
# -----------------------------------------------------------------
#  Example 3: A collection of 2D algebraic benchmark problems
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample3Ext.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
# Benchmark Functions
# -----------------------------------------------------------------
def Ackley(x):
    x = np.atleast_2d(x)
    n = x.shape[1]
    sum_sq = np.sum(x**2, axis=1)
    sum_cos = np.sum(np.cos(2 * np.pi * x), axis=1)
    return -20 * np.exp(-0.2 * np.sqrt(sum_sq / n)) - np.exp(sum_cos / n) + 20 + np.exp(1)

def Beale(x):
    x = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return (1.5 - x1 + x1 * x2)**2 + (2.25 - x1 + x1 * x2**2)**2 + (2.625 - x1 + x1 * x2**3)**2

def Booth(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return (x1 + 2 * x2 - 7)**2 + (2 * x1 + x2 - 5)**2

def BukinN6(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return 100 * np.sqrt(np.abs(x2 - 0.01*x1**2)) + 0.01*np.abs(x1 + 10)

def CrossInTray(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return -0.0001 * (np.abs(np.sin(x1)*np.sin(x2)*np.exp(np.abs(100 - np.sqrt(x1**2+x2**2)/np.pi))) + 1)**0.1

def DixonPrice(x):
    x  = np.atleast_2d(x)
    term1 = (x[:, 0] - 1)**2
    if x.shape[1] > 1:
        sum_seq = np.arange(2, x.shape[1] + 1)
        terms = 2 * sum_seq * ((2 * x[:, 1:]**2 - x[:, :-1])**2)
        return term1 + np.sum(terms, axis=1)
    else:
        return term1

def Easom(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return -np.cos(x1)*np.cos(x2)*np.exp(-((x1-np.pi)**2+(x2-np.pi)**2))

def Eggholder(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return -(x2+47)*np.sin(np.sqrt(np.abs(x2+x1/2+47))) - x1*np.sin(np.sqrt(np.abs(x1-(x2+47))))

def GoldsteinPrice(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    term1 = 1 + (x1 + x2 + 1)**2*(19 - 14*x1 + 3*x1**2 - 14*x2 + 6*x1*x2 + 3*x2**2)
    term2 = 30 + (2*x1 - 3*x2)**2*(18 - 32*x1 + 12*x1**2 + 48*x2 - 36*x1*x2 + 27*x2**2)
    return term1 * term2

def Griewank(x):
    x  = np.atleast_2d(x)
    m, n = x.shape
    i = np.arange(1, n + 1)
    sum_sq = np.sum(x**2, axis=1)/4000
    prod_cos = np.prod(np.cos(x/np.sqrt(i)), axis=1)
    return sum_sq - prod_cos + 1

def Himmelblau(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2

def HolderTable(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return -np.abs(np.sin(x1)*np.cos(x2)*np.exp(np.abs(1 - np.sqrt(x1**2+x2**2)/np.pi)))

def LeviN13(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return np.sin(3*np.pi*x1)**2 + (x1-1)**2*(1+np.sin(3*np.pi*x2)**2) + (x2-1)**2*(1+np.sin(2*np.pi*x2)**2)

def Matyas(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return 0.26*(x1**2+x2**2)-0.48*x1*x2

def McCormick(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return np.sin(x1+x2) + (x1-x2)**2 - 1.5*x1 + 2.5*x2 + 1

def Rastrigin(x):
    x  = np.atleast_2d(x)
    n = x.shape[1]
    A = 10
    return A*n + np.sum(x**2 - A*np.cos(2*np.pi*x), axis=1)

def Rosenbrock(x):
    x  = np.atleast_2d(x)
    return np.sum(100*(x[:, 1:]-x[:, :-1]**2)**2 + (1-x[:, :-1])**2, axis=1)

def SchafferN2(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return 0.5 + (np.sin(x1**2-x2**2)**2 - 0.5)/(1+0.001*(x1**2+x2**2))**2

def SchafferN4(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return 0.5 + (np.cos(np.sin(np.abs(x1**2-x2**2)))**2 - 0.5)/(1+0.001*(x1**2+x2**2))**2

def Shekel(x):
    x  = np.atleast_2d(x)
    m = 10
    A = np.array([[4, 4],
                  [1, 1],
                  [8, 8],
                  [6, 6],
                  [3, 3],
                  [2, 2],
                  [5, 5],
                  [8, 8],
                  [6, 6],
                  [7, 7]])
    C = (1/m)*np.array([1,2,2,4,4,6,3,7,5,5])
    F_val = np.zeros(x.shape[0])
    for i in range(m):
        F_val -= 1.0/(np.sum((x - A[i, :])**2, axis=1) + C[i])
    return F_val

def Sphere(x):
    x  = np.atleast_2d(x)
    return np.sum(x**2, axis=1)

def StyblinskiTang(x):
    x  = np.atleast_2d(x)
    return np.sum(x**4 - 16*x**2 + 5*x, axis=1)/2.0

def ThreeHumpCamel(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    return 2*x1**2 - 1.05*x1**4 + (x1**6)/6 + x1*x2 + x2**2

def Zakharov(x):
    x  = np.atleast_2d(x)
    n = x.shape[1]
    sum1 = np.sum(x**2, axis=1)
    i = np.arange(1, n+1)
    sum2 = np.sum(0.5*i*x, axis=1)
    return sum1 + sum2**2 + sum2**4
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# printExitFlagMeaning - prints interpretation of the exit flag.
# -----------------------------------------------------------------
def printExitFlagMeaning(ExitFlag):

    if ExitFlag == 1:
        print(' Termination: Maximum number of iterations reached.')
    elif ExitFlag == 2:
        print(' Termination: Solution stalled.')
    elif ExitFlag == 3:
        print(' Termination: Maximum number of function evaluations reached.')
    elif ExitFlag == 4:
        print(' Termination: Function value variation is small.')
    elif ExitFlag == 5:
        print(' Termination: Standard deviation variation is small.')
    elif ExitFlag == 6:
        print(' Termination: Minimum function value criterion met.')
    else:
        print(' Termination: Ended with an unknown exit flag.')
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# List of benchmark tests: each entry is a tuple with:
# (function name, function handle, lb, ub, global optima points, global optimum value)
# -----------------------------------------------------------------
benchmarks = [
    ('Ackley', Ackley, [-10, -10], [10, 10], [ [0, 0] ], 0.0),
    ('Beale', Beale, [-5, -5], [5, 5], [ [3, 0.5] ], 0.0),
    ('Booth', Booth, [-20, -20], [20, 20], [ [1, 3] ], 0.0),
    ('BukinN6', BukinN6, [-15, -3], [-5, 3], [ [-10, 1] ], 0.0),
    ('CrossInTray', CrossInTray, [-6, -6], [6, 6],
     [ [-1.3491, -1.3491], [1.3491, 1.3491],
       [-1.3491, 1.3491], [1.3491, -1.3491] ], -2.06261),
    ('DixonPrice', DixonPrice, [-10, -10], [10, 10], [ [1, 1] ], 0.0),
    ('Easom', Easom, [2, 2], [4, 4], [ [3.14159, 3.14159] ], -1.0),
    ('Eggholder', Eggholder, [-512, -512], [512, 512], [ [512, 404.23190] ], -959.64070),
    ('GoldsteinPrice', GoldsteinPrice, [-2, -2], [2, 2], [ [0, -1] ], 3.0),
    ('Griewank', Griewank, [-30, -30], [30, 30], [ [0, 0] ], 0.0),
    ('Himmelblau', Himmelblau, [-5, -5], [5, 5],
     [ [3, 2], [-2.805118, 3.283186],
       [-3.779310, -3.283186], [3.584428, -1.848126] ], 0.0),
    ('HolderTable', HolderTable, [-10, -10], [10, 10],
     [ [8.05502, 9.66459], [-8.05502, 9.66459],
       [8.05502, -9.66459], [-8.05502, -9.66459] ], -19.20850),
    ('LeviN13', LeviN13, [-10, -10], [10, 10], [ [1, 1] ], 0.0),
    ('Matyas', Matyas, [-100, -100], [100, 100], [ [0, 0] ], 0.0),
    ('McCormick', McCormick, [-1.5, -3], [4, 4], [ [-0.54719, -1.54719] ], -1.91330),
    ('Rastrigin', Rastrigin, [-5.12, -5.12], [5.12, 5.12], [ [0, 0] ], 0.0),
    ('Rosenbrock', Rosenbrock, [-2, -1], [2, 3], [ [1, 1] ], 0.0),
    ('SchafferN2', SchafferN2, [-25, -25], [25, 25], [ [0, 0] ], 0.0),
    ('SchafferN4', SchafferN4, [-25, -25], [25, 25],
     [ [0, 1.25313], [0, -1.25313], [1.25313, 0], [-1.25313, 0] ], 0.29258),
    ('Shekel', Shekel, [-2, -2], [12, 12], [ [4, 4] ], -10.15320),
    ('Sphere', Sphere, [-5, -5], [5, 5], [ [0, 0] ], 0.0),
    ('StyblinskiTang', StyblinskiTang, [-5, -5], [5, 5], [ [-2.90353, -2.90353] ], -78.33230),
    ('ThreeHumpCamel', ThreeHumpCamel, [-5, -5], [5, 5], [ [0, 0] ], 0.0),
    ('Zakharov', Zakharov, [-5, -5], [10, 10], [ [0, 0] ], 0.0)
]
# -----------------------------------------------------------------

# Custom colormap and colors
Ncolors     = 20
MyColorMap1 = colormaps['viridis'].resampled(Ncolors)

MyYellow    = [0.9290, 0.6940, 0.1250]
MyBlue      = [0.0000, 0.4470, 0.7410]
MyRed       = [0.6350, 0.0780, 0.1840]

# Loop over each benchmark test
for bench in benchmarks:
    # CEopt customization
    CEstr = {}
    CEstr['Verbose'] = False
    CEstr['Nsamp'  ] = 1000

    # Force vectorized function calls to avoid 1D inputs.
    CEstr['isVectorized'] = False
    
    # Extract benchmark info
    funcName, funcHandle, lb_bench, ub_bench, globalOptimaPoints, globalOptimaValue = bench
    lb_bench = np.array(lb_bench)
    ub_bench = np.array(ub_bench)
    Nvars    = lb_bench.size
    
    # Initialize mean and standard deviation vectors
    mu0 = lb_bench + (ub_bench - lb_bench) * np.random.rand(Nvars)
    sigma0 = (ub_bench - lb_bench) / 2.0
    
    # Define tolerance level for convergence
    tol = 0.5
    
    print('\n---------------------------------------------------')
    print(' Optimizing ' + funcName + ' Function...')
    
    start_time = time.time()
    Xopt, Fopt, ExitFlag, CEstr = CEopt(funcHandle, mu0, sigma0, lb_bench, ub_bench, None, CEstr)
    elapsedTime = time.time() - start_time
    
    # Find the closest analytical optimum
    minDist = np.inf
    closestOptimum = None
    for opt in globalOptimaPoints:
        opt_arr = np.array(opt)
        dist    = np.linalg.norm(Xopt - opt_arr)
        if dist < minDist:
            minDist        = dist
            closestOptimum = opt_arr
    isCloseToOptimum = minDist < tol
    
    print(' Numerical Solution:')
    print('  Xopt: ' + np.array2string(Xopt, precision=5))
    print('  Fopt: ' + str(Fopt))
    print(' Analytical Solution:')
    for opt in globalOptimaPoints:
        print('  X*: ' + np.array2string(np.array(opt), precision=5) + ', F*: ' + str(globalOptimaValue))
    print(' Closest global optimum: ' + np.array2string(closestOptimum, precision=5))
    print(' Minimum distance to global optimum: ' + str(minDist))
    printExitFlagMeaning(ExitFlag)
    if isCloseToOptimum:
        print(' CEopt successfully approximated a global optimum.')
    else:
        print(' CEopt did not approximate a global optimum within the specified tolerance.')
    print('---------------------------------------------------')
    print('Elapsed time: {:.4f} seconds'.format(elapsedTime))
    
    # Generate meshgrid for contour plot
    Ngrid  = 100
    xRange = np.linspace(lb_bench[0], ub_bench[0], Ngrid)
    yRange = np.linspace(lb_bench[1], ub_bench[1], Ngrid)
    X, Y   = np.meshgrid(xRange, yRange)
    pts    = np.column_stack((X.ravel(), Y.ravel()))
    Z      = funcHandle(pts)
    Z      = Z.reshape(X.shape)
    
    plt.figure(figsize=(8, 6))
    plt.clf()
    plt.contourf(X, Y, Z, 30)
    plt.set_cmap(MyColorMap1)
    plt.colorbar()
    plt.title(funcName + ' Function')
    plt.xlabel('x1', fontsize=20, fontname='Arial')
    plt.ylabel('x2', fontsize=20, fontname='Arial')
    
    # Plot analytical global optima
    for idx, opt in enumerate(globalOptimaPoints):
        opt_arr = np.array(opt)
        if idx == 0:
            plt.plot(opt_arr[0], opt_arr[1], 'o', color=MyYellow, markersize=15, label='Analytical Solution')
        else:
            plt.plot(opt_arr[0], opt_arr[1], 'o', color=MyYellow, markersize=15, label='_nolegend_')
    
    # Plot initial approximation
    plt.plot(mu0[0], mu0[1], 'd', color=MyRed, linewidth=1, markersize=7,
             markerfacecolor=MyRed, label='Initial Approximation')
    
    # Plot CE iterations
    for j in range(CEstr['iter']):
        if j == 0:
            plt.plot(CEstr['xmean'][j, 0], CEstr['xmean'][j, 1], 'd', color=MyRed,
                     linewidth=1, markersize=7, label='CE Iterations')
        else:
            plt.plot(CEstr['xmean'][j, 0], CEstr['xmean'][j, 1], 'd', color=MyRed,
                     linewidth=1, markersize=7, label='_nolegend_')
    
    # Plot CE solution
    plt.plot(Xopt[0,0], Xopt[0,1], 'x', color=MyBlue, linewidth=3, markersize=25, label='CE Optimum')
    
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontname('Arial')
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    ax.tick_params(labelsize=18)
    
    # Use a valid legend location ('lower left' instead of 'southwest')
    plt.legend(loc='lower left', fontsize=12)
    plt.box(True)
    
    # Adjust subplot margins manually to help tight_layout succeed
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
    
    saveFileName = 'ContourPlot_' + funcName + '.png'
    plt.savefig(saveFileName, format='png', dpi=300)

    # Close the current figure
    plt.close()
