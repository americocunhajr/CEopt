# -----------------------------------------------------------------
#  CEopt.py
# -----------------------------------------------------------------
#  Programmer: Americo Cunha Jr
#              americo.cunhajr@gmail.com
# 
#  Originally programmed in: Mar 27, 2025
#           Last updated in: Jun 12, 2025
# -----------------------------------------------------------------
#  This routine employs the Cross-entropy (CE) method to solve the 
#  following optimization problem:
# 
#                       Xopt = arg min F(x)
#  
#                       subject to
# 
#                       lb <= x <= ub
#                       g(x) <= 0
#                       h(x)  = 0
#  where:
#  - F   : R^Nvars -> R is a given scalar objective function
#  - lb  : (1 x Nvars) vector of lower bounds for the decision variables
#  - ub  : (1 x Nvars) vector of upper bounds for the decision variables
#  - g(x): R^Nvars -> R^m is a vector of m inequality constraint functions
#  - h(x): R^Nvars -> R^p is a vector of p equality   constraint functions
#
#  The goal is to minimize the scalar objective function F(x) defined 
#  within a known rectangular domain (feasible region), while also 
#  satisfying the given equality and inequality constraints. The feasible 
#  region is further defined by the constraints g(x) <= 0 and h(x) = 0, 
#  in addition to the bounds lb <= x <= ub.
#
#  The algorithm samples the feasible region using a truncated Gaussian
#  distribution and updates its parameters with the aid of an elite set
#  (defined by the better samples), aiming to transform this Gaussian
#  distribution into a Dirac distribution centered at the global optimum.
#  For constrained optimization problems, the CE method integrates a 
#  mechanism for handling the constraints, such as the augmented 
#  Lagrangian method, to incorporate the effects of constraints into the 
#  optimization process.
#
#  Input:
#  fun    - Function handle for the objective function. This function
#           must accept a 1 x Nvars row vector (representing a single 
#           sample) or an M x Nvars matrix (representing M samples with 
#           variables in columns) as input and return a scalar value or 
#           a row vector of M scalar values (for vectorized operations) 
#           respectively.
#  xmean0  - (1 x Nvars) initial mean
#  sigma0  - (1 x Nvars) initial standard deviation
#  lb      - (1 x Nvars) lower bound
#  ub      - (1 x Nvars) upper bound
#  nonlcon - Function handle for thenonlinear constraint function.
#  CEstr   - Struct (here, a dictionary) containing parameters and settings for the CE method.
# 
#  CEstr fields include:
#  * Verbose          : boolean flag to enable/disable screen output
#  * isConstrained    : boolean flag to indicate a constrained problem
#  * isVectorized     : boolean flag to indicate a vectorized function
#  * Nvars            : number of design variables in x
#  * EliteFactor      : proportion of samples for the elite set
#  * Nsamp            : number of samples to draw per iteration
#  * MaxIter          : maximum number of iterations
#  * MaxStall         : maximum number of stall iterations
#  * MaxFcount        : maximum number of objective function evaluations
#  * MinFval          : minimum admissible value for objective function 
#  * TolAbs           : absolute tolerance
#  * TolRel           : relative tolerance
#  * TolCon           : constraint violation tolerance
#  * TolFun           : function value change tolerance
#  * alpha            : smoothing parameter for the mean update
#  * beta             : smoothing parameter for the std. dev. update
#  * q                : exponent in dynamic smoothing parameter
#  * NonlconAlgorithm : algorithm to handle nonlinear constraints
#  * InitialPenalty   : initial penalty value for the augmented Lagrangian
#  * PenaltyFactor    : factor by which the penalty parameter is increased
#  * MaximumPenalty   : maximum value for the penalty parameter
#  * xmean            : history of mean value over iterations
#  * xmedian          : history of median over iterations
#  * xbest            : history of best sample point over iterations
#  * Fmean            : history of objective function mean over elite set
#  * Fmedian          : history of objective function median over elite set
#  * Fbest            : history of objective function best value found
#  * sigma            : history of standard deviation over iterations
#  * ErrorS           : history of standard deviation error
#  * ErrorC           : history of constraint violation error (for constrained problems)
#  * iter             : total number of iterations performed
#  * stall            : number of iterations without significant progress
#  * Fcount           : total number of function evaluations performed
#  * ConvergenceStatus: boolean flag indicating if the algorithm converged
#                       (true) or not (false).
#            
#  This struct can also include additional fields for customized behavior 
#  or extensions of the CE method.
# 
#  Output:
#  Xopt     - (1 x Nvars) optimal point
#  Fopt     - scalar optimal value
#  ExitFlag - Flag indicating the reason for algorithm termination:
#             0 - algorithm is running or has not been initialized
#             1 - maximum number of iterations reached
#             2 - no significant change in objective function (stalled)
#             3 - maximum number of function evaluations reached
#             4 - function change and constraint error smaller than tolerance
#             5 - std dev and constraint errors smaller than tolerance
#             6 - minimum function value criterion met
#  CEstr   - The updated Cross-Entropy object struct containing the final
#            state of the algorithm and possibly additional diagnostic information.
# -----------------------------------------------------------------
#  References:
# 
#  [1] Reuven Y. Rubinstein, Dirk P. Kroese,
#      The Cross-Entropy Method: A Unified Approach to Combinatorial 
#      Optimization, Monte-Carlo Simulation, and Machine Learning,
#      Springer-Verlag, 2004.
# 
#  [2] A. Cunha Jr, M. V. Issa, J. C. Basilio, J. G. Telles Ribeiro,
#      CEopt: A MATLAB Package for Non-convex Optimization with the
#      Cross-Entropy Method, ArXiv, 2024
# -----------------------------------------------------------------
#  Copyright (C) 2025  Americo Cunha Jr et al.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# Imported libraries
# -----------------------------------------------------------------
import numpy as np
from scipy.special import erfc, erfcinv
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# CEopt - Cross-Entropy optimization method main function
# -----------------------------------------------------------------
def CEopt(fun, xmean0, sigma0, lb, ub, nonlcon=None, CEstr=None):

    # consistency check for the mandatory inputs parameters
    lb,ub,xmean0,sigma0,Nvars = CheckInput(lb,ub,xmean0,sigma0)
    
    # check if CEstr is not provided or is empty
    if CEstr is None or not CEstr:
        CEstr = {}
    
    # set the default parameters for CEstr (if necessary)
    CEstr = InitializeCEstr(CEstr,Nvars)
    
    # consistency check for CEstr parameters
    CheckCEstr(CEstr)
    
    # check if nonlcon is not provided or is empty
    if nonlcon is None:
        CEstr['isConstrained'] = False
        # objective function
        ObjFun = lambda x, *args: fun(x)
    else:
        CEstr['isConstrained'] = True
        # objective function
        if CEstr['NonlconAlgorithm'] == 'AugLagLog':
            ObjFun = lambda x, Lg, Lh, p: AugLagrangian1(x,fun,nonlcon,Lg,Lh,p)
        else:
            ObjFun = lambda x, Lg, Lh, p: AugLagrangian2(x,fun,nonlcon,Lg,Lh,p)
    
    # decide the appropriate optimization solver
    if not CEstr['isConstrained']:
        # CE solver for an unconstrained problem
        Xopt,Fopt,ExitFlag,CEstr = UnconstrSolverCE(ObjFun,Nvars,xmean0,sigma0,lb,ub,CEstr)
    else:
        # CE solver for a constrained problem
        Xopt,Fopt,ExitFlag,CEstr = ConstrSolverCE(ObjFun,Nvars,xmean0,sigma0,lb,ub,nonlcon,CEstr)
    
    return Xopt, Fopt, ExitFlag, CEstr
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# CheckInput - verify input parameters for possible errors
# -----------------------------------------------------------------
def CheckInput(lb, ub, xmean0, sigma0):

    # check if lb and ub are empty
    if lb is None or ub is None:
        raise ValueError('lb and ub must be non-empty')

    # ensure lb and ub are row vectors
    lb = np.array(lb, ndmin=1).reshape(1, -1)
    ub = np.array(ub, ndmin=1).reshape(1, -1)

    # check if lb and ub are numpy arrays
    if not isinstance(lb, np.ndarray):
        raise TypeError('lb must be a numpy array')
    if not isinstance(ub, np.ndarray):
        raise TypeError('ub must be a numpy array')
    
    # check if lb and ub are empty vectors
    if lb.size == 0 or ub.size == 0:
        raise ValueError('lb and ub must be non-empty vectors')
    
    # number of variables
    Nvars = int(lb.size)
    
    # check for consistency in lb and ub
    if ub.size != Nvars:
        raise ValueError('lb and ub must have the same dimension')
    if np.isnan(lb).any() or np.isnan(ub).any():
        raise ValueError('lb and ub cannot have a NaN components')
    if np.any(lb >= ub):
        raise ValueError('lb < ub for all components')
    
    # define xmean0 (if necessary)
    if xmean0 is None or (isinstance(xmean0, np.ndarray) and xmean0.size == 0):
        xmean0 = (ub + lb) / 2.0
    else:
        xmean0 = np.array(xmean0, ndmin=1).reshape(1, -1)
        if not isinstance(xmean0, np.ndarray):
            raise TypeError('xmean0 must be a numpy array')
        if xmean0.size == 0:
            xmean0 = (ub + lb) / 2.0
    
    # define sigma0 (if necessary)
    if sigma0 is None or (isinstance(sigma0, np.ndarray) and sigma0.size == 0):
        sigma0 = (ub - lb) / np.sqrt(12.0)
    else:
        sigma0 = np.array(sigma0, ndmin=1).reshape(1, -1)
        if not isinstance(sigma0, np.ndarray):
            raise TypeError('sigma0 must be a numpy array')
        if sigma0.size == 0:
            sigma0 = (ub - lb) / np.sqrt(12.0)
    
    # check for consistency in xmean0 and sigma0
    if xmean0.shape != sigma0.shape:
        raise ValueError('xmean0 and sigma0 must have the same dimensions')
    if np.isnan(xmean0).any() or np.isnan(sigma0).any():
        raise ValueError('xmean0 and sigma0 cannot have a NaN components')
    if np.isinf(xmean0).any() or np.isinf(sigma0).any():
        raise ValueError('xmean0 and sigma0 cannot have an Inf components')
    if np.any(xmean0 < lb) or np.any(xmean0 > ub):
        raise ValueError('xmean0 must be in [lb,ub] interval')
    if np.any(sigma0 <= 0.0):
        raise ValueError('All components of sigma0 must be positive')
    
    # check for dimension consistency in xmean0 and sigma0
    if xmean0.size != Nvars:
        raise ValueError('xmean0 must be a 1 x Nvars vector')
    if sigma0.size != Nvars:
        raise ValueError('sigma0 must be a 1 x Nvars vector')
    
    return lb, ub, xmean0, sigma0, Nvars
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# InitializeCEstr - initialize and set default parameters for CEstr
# -----------------------------------------------------------------
def InitializeCEstr(CEstr, Nvars):
    DefaultParams = {
        'Verbose'         : True,
        'isConstrained'   : False,
        'isVectorized'    : False,
        'Nvars'           : Nvars,
        'EliteFactor'     : 0.05,
        'Nsamp'           : 100,
        'MaxIter'         : 100 * Nvars,
        'MaxStall'        : 50,
        'MaxFcount'       :  np.inf,
        'MinFval'         : -np.inf,
        'TolAbs'          : 1.0e-6,
        'TolRel'          : 1.0e-3,
        'TolCon'          : 1.0e-3,
        'TolFun'          : 1.0e-3,
        'alpha'           : 0.4,
        'beta'            : 0.4,
        'q'               : 10.0,
        'NonlconAlgorithm': 'AugLagLog',
        'InitialPenalty'  : 10.0,
        'PenaltyFactor'   : 10.0,
        'MaximumPenalty'  : np.inf
    }
    
    # assign values for undefined fields in CEstr
    for key, val in DefaultParams.items():
        if key not in CEstr:
            CEstr[key] = val

    # Remove strange fields
    allowed_keys = set(DefaultParams.keys())
    CEstr = {k: v for k, v in CEstr.items() if k in allowed_keys}
    
    # # get the strange fields in CEstr (fields not in DefaultParams)
    # StrangeFields = {k: v for k, v in CEstr.items() if k not in DefaultParams}

    # # remove the strange fields from CEstr
    # for k in list(StrangeFields.keys()):
    #     CEstr.pop(k)

    # order the default fields in CEstr
    CEstr_ordered = {k: CEstr[k] for k in DefaultParams.keys() if k in CEstr}
    CEstr = CEstr_ordered

    # preallocate memory for histories (using np.full with np.nan)
    maxiter          = CEstr['MaxIter']
    CEstr['xmean'  ] = np.full((maxiter, Nvars), np.nan)
    CEstr['xmedian'] = np.full((maxiter, Nvars), np.nan)
    CEstr['xbest'  ] = np.full((maxiter, Nvars), np.nan)
    CEstr['Fmean'  ] = np.full((maxiter, 1    ), np.nan)
    CEstr['Fmedian'] = np.full((maxiter, 1    ), np.nan)
    CEstr['Fbest'  ] = np.full((maxiter, 1    ), np.nan)
    CEstr['sigma'  ] = np.full((maxiter, Nvars), np.nan)
    CEstr['ErrorS' ] = np.full((maxiter, 1    ), np.nan)
    CEstr['ErrorC' ] = np.full((maxiter, 1    ), np.nan)
    CEstr['iter'   ] = 0
    CEstr['stall'  ] = 0
    CEstr['Fcount' ] = 0
    
    return CEstr
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# is_number - check if the input is of float type
# -----------------------------------------------------------------
def is_number(x):
    return isinstance(x, (int, float, np.integer, np.floating))

# -----------------------------------------------------------------
# is_integer - check if the input is of integer type
# -----------------------------------------------------------------
def is_integer(x):
    return isinstance(x, (int, np.integer)) or (isinstance(x, float) and x.is_integer())

# -----------------------------------------------------------------
# CheckCEstr - check parameters consistency for CEstr
# -----------------------------------------------------------------
def CheckCEstr(CEstr):

    # Verbose must be boolean
    if not isinstance(CEstr['Verbose'], bool):
        raise ValueError('Verbose must be boolean')

    # EliteFactor check
    if (not is_number(CEstr['EliteFactor']) or CEstr['EliteFactor'] <= 0 or CEstr['EliteFactor'] > 1):
        raise ValueError('EliteFactor must be such that 0 < EliteFactor <= 1')

    # Nsamp check
    if not is_number(CEstr['Nsamp']) or CEstr['Nsamp'] <= 1:
        raise ValueError('Nsamp must be greater than 1')

    # MaxIter check
    if not is_integer(CEstr['MaxIter']) or CEstr['MaxIter'] < 1:
        raise ValueError('MaxIter must be a positive integer')

    # MaxStall check
    if not is_integer(CEstr['MaxStall']) or CEstr['MaxStall'] < 1:
        raise ValueError('MaxStall must be a positive integer')

    # MaxFcount check
    fc = CEstr['MaxFcount']
    if (is_integer(fc) and fc < 1) or (not is_integer(fc) and fc != np.inf):
        raise ValueError('MaxFcount must be a positive integer or infinity')

    # MinFval check
    if not is_number(CEstr['MinFval']):
        raise ValueError('MinFval must be numeric')

    # TolAbs check
    if not is_number(CEstr['TolAbs']) or CEstr['TolAbs'] <= 0.0:
        raise ValueError('TolAbs must be positive real')

    # TolRel check
    if not is_number(CEstr['TolRel']) or CEstr['TolRel'] < 0.0:
        raise ValueError('TolRel must be non-negative real')

    # TolCon check
    if not is_number(CEstr['TolCon']) or CEstr['TolCon'] < 0.0:
        raise ValueError('TolCon must be non-negative real')

    # TolFun check
    if not is_number(CEstr['TolFun']) or CEstr['TolFun'] < 0.0:
        raise ValueError('TolFun must be non-negative')

    # alpha check
    if not is_number(CEstr['alpha']) or CEstr['alpha'] <= 0 or CEstr['alpha'] > 1:
        raise ValueError('alpha must be such that 0 < alpha <= 1')

    # beta check
    if not is_number(CEstr['beta']) or CEstr['beta'] <= 0:
        raise ValueError('beta must be non-negative')

    # q check
    if not is_number(CEstr['q']) or CEstr['q'] <= 0:
        raise ValueError('q must be non-negative')

    # NonlconAlgorithm check
    if not isinstance(CEstr['NonlconAlgorithm'], str) or \
       CEstr['NonlconAlgorithm'] not in ['AugLagLog', 'AugLagMax']:
        raise ValueError('Unknown option for NonlconAlgorithm')

    # InitialPenalty check
    if not is_number(CEstr['InitialPenalty']) or CEstr['InitialPenalty'] <= 0:
        raise ValueError('InitialPenalty must be non-negative')

    # PenaltyFactor check
    if not is_number(CEstr['PenaltyFactor']) or CEstr['PenaltyFactor'] <= 1:
        raise ValueError('PenaltyFactor must be greater than 1')

    if not is_number(CEstr['MaximumPenalty']) or (CEstr['MaximumPenalty'] <= 1 and CEstr['MaximumPenalty'] != np.inf):
        raise ValueError('MaximumPenalty must be greater than 1 or infinity')
# -----------------------------------------------------------------
    

# -----------------------------------------------------------------
# UnconstrSolverCE - solve an unconstrained optimization problem
# -----------------------------------------------------------------
def UnconstrSolverCE(fun, Nvars, xmean0, sigma0, lb, ub, CEstr):
    
    t           = 0                                # iteration counter
    stall       = 0                                # stall iterations counter
    Fcount      = 0                                # function evaluation counter
    EliteFactor = CEstr['EliteFactor']             # elite factor
    Nsamp       = CEstr['Nsamp'      ]             # total number of samples
    Nelite      = int(round(EliteFactor*Nsamp))    # number of elite samples
    MaxIter     = CEstr['MaxIter'    ]             # maximum iterations
    TolAbs      = CEstr['TolAbs'     ]             # absolute tolerance
    TolRel      = CEstr['TolRel'     ]             # relative tolerance
    alpha       = CEstr['alpha'      ]             # smoothing parameter for mean
    beta        = CEstr['beta'       ]             # smoothing parameter for std. dev.
    q           = CEstr['q'          ]             # dynamic update parameter
    Xopt        = np.full_like(xmean0, np.nan)     # optimal point
    Fopt        = np.inf                           # optimal value
    ExitFlag    = 0                                # termination condition flag

    # preallocate memory for design variables samples
    X = np.zeros((Nsamp, Nvars))
    
    # preallocate memory for objective function evaluations
    F = np.empty((Nsamp, 1))*np.nan

    # one-time dimensions validation flag
    FirstCheckFlag = True
    
    # loop to sample the domain and update the distribution
    while ExitFlag == 0 and t <= MaxIter:

        # update level counter
        t += 1
        
        # sample the domain from a truncated Gaussian distribution
        X = DomainSampling(xmean0,sigma0,lb,ub,Nvars,Nsamp,X)

        # evaluate objective function
        if not CEstr['isVectorized']:
            # case where fun is not a vectorized function
            for n in range(Nsamp):
                F[n,0] = fun(X[n, :].reshape(1,-1))
        else:
            # case where fun is a vectorized function
            F = fun(X)
            if FirstCheckFlag:
                if not isinstance(F,np.ndarray) or F.shape != (Nsamp,1):
                    raise ValueError('Vectorized function must return Nsamp x 1 array')
                FirstCheckFlag = False
        
        # update function evaluation counter
        Fcount += Nsamp
        
        # define elite samples set 
        EliteSetId = DefineEliteSet(F, Nelite)
        
        # update the distribution parameters
        xmean,xmedian,xbest,Fmean,Fmedian,Fbest,sigma = UpdateDistribution(F,X,EliteSetId,xmean0,sigma0,alpha,beta,q,t)
        
        # update standard deviation error
        ErrorS, SmallErrorS = ComputeErrorS(sigma, sigma0, TolAbs, TolRel)
        
        # update old parameters
        xmean0 = xmean.copy()
        sigma0 = sigma.copy()
        
        # update the optimum
        if Fbest < Fopt:
            CEstr['xbest'][t - 1, :] = xbest.copy()
            CEstr['Fbest'][t - 1, 0] = Fbest
            Xopt                     = xbest.copy()
            Fopt                     = Fbest
            stall                    = 0
        else:
            CEstr['xbest'][t - 1, :] = CEstr['xbest'][t - 2, :].copy()
            CEstr['Fbest'][t - 1, 0] = CEstr['Fbest'][t - 2, 0]
            stall                   += 1
        
        # update optimization process history
        CEstr['iter'   ]           = t
        CEstr['stall'  ]           = stall
        CEstr['Fcount' ]           = Fcount
        CEstr['xmean'  ][t - 1, :] = xmean.copy()
        CEstr['xmedian'][t - 1, :] = xmedian.copy()
        CEstr['Fmean'  ][t - 1, 0] = Fmean
        CEstr['Fmedian'][t - 1, 0] = Fmedian
        CEstr['sigma'  ][t - 1, :] = sigma.copy()
        CEstr['ErrorS' ][t - 1, 0] = ErrorS
        
        # print iteration progress on the screen
        if CEstr['Verbose']:
            PrintProgress(t, Nvars, CEstr)
        
        # check the convergence
        ExitFlag = CheckConv(Fopt,SmallErrorS,None,CEstr)
    
    # convergence check and update of 'ConvergenceStatus' field
    if ExitFlag > 3:
        CEstr['ConvergenceStatus'] = True
    else:
        CEstr['ConvergenceStatus'] = False
    
    # print resume
    if CEstr['Verbose']:
        PrintEnd(Xopt, Fopt, ExitFlag, CEstr)
    
    # delete empty entries from sampling records
    CEstr = DeleteEmptyEntries(t, CEstr)
    
    return Xopt, Fopt, ExitFlag, CEstr
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# ConstrSolverCE - solve a constrained optimization problem
# -----------------------------------------------------------------
def ConstrSolverCE(fun, Nvars, xmean0, sigma0, lb, ub, nonlcon, CEstr):

    t           = 0                             # iteration counter
    stall       = 0                             # stall iterations counter
    Fcount      = 0                             # function evaluation counter
    EliteFactor = CEstr['EliteFactor']          # elite factor
    Nsamp       = CEstr['Nsamp'      ]          # total number of samples
    Nelite      = int(round(EliteFactor*Nsamp)) # number of elite samples
    MaxIter     = CEstr['MaxIter'    ]          # maximum iterations
    TolAbs      = CEstr['TolAbs'     ]          # absolute tolerance
    TolRel      = CEstr['TolRel'     ]          # relative tolerance
    TolCon      = CEstr['TolCon'     ]          # constraint tolerance
    alpha       = CEstr['alpha'      ]          # smoothing parameter for mean
    beta        = CEstr['beta'       ]          # smoothing parameter for std. dev.
    q           = CEstr['q'          ]          # dynamic update parameter
    Xopt        = np.full_like(xmean0,np.nan)   # optimal point
    Fopt        = np.inf                        # optimal value
    ExitFlag    = 0                             # termination condition flag
    
    # initialize penalty parameters
    Penalty        = CEstr['InitialPenalty']
    PenaltyFactor  = CEstr['PenaltyFactor' ]
    MaximumPenalty = CEstr['MaximumPenalty']
    
    # initialize Lagrange multipliers
    #   G0      - (1 x Ng) inequality constraints vector at xmean0
    #   H0      - (1 x Nh) equality   constraints vector at xmean0
    #   lambdaG - (1 x Ng) Lagrange multipliers for inequality constraints
    #   lambdaH - (1 x Nh) Lagrange multipliers for   equality constraints
    G0, H0 = nonlcon(xmean0)              # xmean0: (1, Nvars)
    if G0 is None or np.size(G0) == 0:
        G0 = np.zeros((1, 1))
    if H0 is None or np.size(H0) == 0:
        H0 = np.zeros((1, 1))
    lambdaG = np.zeros_like(G0)
    lambdaH = np.zeros_like(H0)
    
    # initialize constraint error
    ErrorC, _ = ComputeErrorC(G0,H0,lambdaG,lambdaH,Penalty,TolCon,1.0)
    
    # preallocate memory for design variables samples
    X = np.zeros((Nsamp, Nvars))
    
    # preallocate memory for objective function evaluations
    AL = np.empty((Nsamp, 1))*np.nan
    F  = np.empty((Nsamp, 1))*np.nan

    # one-time dimensions validation flag
    FirstCheckFlag = True
    
    # loop to sample the domain and update the distribution
    while ExitFlag == 0 and t <= MaxIter:

        # update level counter
        t += 1
        
        # sample the domain from a truncated Gaussian distribution
        X = DomainSampling(xmean0,sigma0,lb,ub,Nvars,Nsamp,X)
        
        # evaluate augmented Lagrangian 
        if not CEstr['isVectorized']:
            # case where fun is not a vectorized function
            for n in range(Nsamp):
                AL[n,0], F[n,0], _, _ = fun(X[n,:].reshape(1,-1),lambdaG,lambdaH,Penalty)
        else:
            # case where fun is a vectorized function
            AL, F, _, _ = fun(X,lambdaG,lambdaH,Penalty)
            if FirstCheckFlag:
                if not isinstance(F,np.ndarray) or F.shape != (Nsamp,1):
                    raise ValueError('Vectorized function must return Nsamp x 1 NumPy array')
                FirstCheckFlag = False
        
        # update function evaluation counter
        Fcount += Nsamp
        
        # define elite samples set 
        EliteSetId = DefineEliteSet(AL, Nelite)
        
        # update the distribution parameters
        xmean,xmedian,xbest,Fmean,Fmedian,Fbest,sigma = UpdateDistribution(AL,X,EliteSetId,xmean0,sigma0,alpha,beta,q,t)
        
        # standard deviation error
        ErrorS,SmallErrorS = ComputeErrorS(sigma,sigma0,TolAbs,TolRel)
        
        # evaluate the constraints at xbest
        G, H = nonlcon(xbest)
        if G is None or np.size(G) == 0:
            G = np.zeros(lambdaG.shape)
        if H is None or np.size(H) == 0:
            H = np.zeros(lambdaH.shape) 
        
        # evaluate the constraint error
        ErrorC,SmallErrorC = ComputeErrorC(G,H,lambdaG,lambdaH,Penalty,TolCon,ErrorC)
        
        # update Lagrange multipliers
        lambdaG,lambdaH = UpdateLagrangeMult(G,H,lambdaG,lambdaH,Penalty)
        
        # update penalty parameter
        Penalty = UpdatePenalty(Penalty,PenaltyFactor,MaximumPenalty,SmallErrorC)
        
        # update old parameters
        xmean0 = xmean.copy()
        sigma0 = sigma.copy()
        
        # update the optimum
        if Fbest < Fopt:
            Fbest                    = F[EliteSetId[:, 0]][0]
            CEstr['xbest'][t - 1, :] = xbest.copy()
            CEstr['Fbest'][t - 1, 0] = Fbest
            Xopt                     = xbest.copy()
            Fopt                     = Fbest
            stall                    = 0
        else:
            CEstr['xbest'][t - 1, :] = CEstr['xbest'][t - 2, :].copy()
            CEstr['Fbest'][t - 1, 0] = CEstr['Fbest'][t - 2, 0]
            stall                   += 1
        
        # update optimization process history
        CEstr['iter'   ]           = t
        CEstr['stall'  ]           = stall
        CEstr['Fcount' ]           = Fcount
        CEstr['xmean'  ][t - 1, :] = xmean.copy()
        CEstr['xmedian'][t - 1, :] = xmedian.copy()
        CEstr['Fmean'  ][t - 1, 0] = Fmean
        CEstr['Fmedian'][t - 1, 0] = Fmedian
        CEstr['sigma'  ][t - 1, :] = sigma.copy()
        CEstr['ErrorS' ][t - 1, 0] = ErrorS
        CEstr['ErrorC' ][t - 1, 0] = ErrorC
        
        # print iteration progress on the screen
        if CEstr['Verbose']:
            PrintProgress(t, Nvars, CEstr)
        
        # check the convergence
        ExitFlag = CheckConv(Fopt,SmallErrorS,SmallErrorC,CEstr)
    
    # convergence check and update of 'ConvergenceStatus' field
    if ExitFlag > 3:
        CEstr['ConvergenceStatus'] = True
    else:
        CEstr['ConvergenceStatus'] = False
    
    # print resume
    if CEstr['Verbose']:
        PrintEnd(Xopt, Fopt, ExitFlag, CEstr)
    
    # delete empty entries from sampling records
    CEstr = DeleteEmptyEntries(t, CEstr)
    
    return Xopt, Fopt, ExitFlag, CEstr
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# AugLagrangian1 - augmented Lagrangian function (log formulation)
#  lambdaG - (1 x Ng) Lagrange multipliers for inequality constraints
#  lambdaH - (1 x Nh) Lagrange multipliers for   equality constraints
# -----------------------------------------------------------------
def AugLagrangian1(x, fun, nonlcon, lambdaG, lambdaH, Penalty):

    # Evaluate objective and constraints
    #   F - (Nsamp x  1) objective function samples
    #   G - (Nsamp x Ng) inequality constraints samples
    #   H - (Nsamp x Nh)   equality constraints samples
    F    = fun(x)
    G, H = nonlcon(x)

    # Handle empty constraints
    if G is None or np.size(G) == 0:
        G = np.zeros((x.shape[0], 1))
    if H is None or np.size(H) == 0:
        H = np.zeros((x.shape[0], 1))

    # Machine epsilon
    eps = np.finfo(float).eps

    # Shift
    s = lambdaG / Penalty

    # Compute the augmented Lagrangian
    term1 = -np.sum(np.multiply(s*lambdaG,np.log(s-G+eps)),axis=1).reshape(-1, 1)
    term2 =  np.dot(H,lambdaH.T)
    term3 = 0.5*Penalty*np.sum(                       H**2,axis=1).reshape(-1, 1)
    #term2 =  np.sum(H*lambdaH                             ,axis=1).reshape(-1, 1)
    #term2 =  np.sum(H*lambdaH, axis=1, keepdims=True)
    AL    = F + term1 + term2 + term3

    return AL, F, G, H
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# AugLagrangian2 - augmented Lagrangian function (max formulation)
#  lambdaG - (1 x Ng) Lagrange multipliers for inequality constraints
#  lambdaH - (1 x Nh) Lagrange multipliers for   equality constraints
# -----------------------------------------------------------------
def AugLagrangian2(x, fun, nonlcon, lambdaG, lambdaH, Penalty):

    # Evaluate objective and constraints
    #   F - (Nsamp x  1) objective function samples
    #   G - (Nsamp x Ng) inequality constraints samples
    #   H - (Nsamp x Nh)   equality constraints samples
    F    = fun(x)
    G, H = nonlcon(x)

    # Handle empty constraints
    if G is None or np.size(G) == 0:
        G = np.zeros((x.shape[0], 1))
    if H is None or np.size(H) == 0:
        H = np.zeros((x.shape[0], 1))

    # Compute shifted constraints.
    H_s = H + lambdaH/Penalty
    G_s = G + lambdaG/Penalty

    # Compute the augmented Lagrangian.
    term1 = np.sum(              H_s**2,axis=1).reshape(-1, 1)
    term2 = np.sum(np.maximum(0,G_s)**2,axis=1).reshape(-1, 1)
    AL    = F + 0.5*Penalty*(term1+term2)

    return AL, F, G, H
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# DefineEliteSet - define elite samples set
# -----------------------------------------------------------------
def DefineEliteSet(F, Nelite):

    # sort objective function evaluations (order statistics)
    Isort = np.argsort(F,axis=0)

    # elite samples indices
    EliteSetId = Isort[:Nelite].reshape(-1, 1)

    return EliteSetId
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# UpdateDistribution - update distribution parameters
# -----------------------------------------------------------------
def UpdateDistribution(F, X, EliteSetId, xmean0, sigma0, alpha, beta, q, t):

    # elite samples and values
    Felite = F[EliteSetId[:, 0]]
    Xelite = X[EliteSetId[:, 0], :]
    
    # estimators for the objective function minimum point
    xmean   = np.mean(Xelite, axis=0).reshape(1,-1)
    xmedian = np.median(Xelite, axis=0).reshape(1,-1)
    xbest   = Xelite[0, :].reshape(1,-1)
    
    # estimators for the objective function minimum value
    Fmean = np.mean(Felite)
    Fmedian = np.median(Felite)
    Fbest = Felite[0]
    
    # estimator for the standard deviation
    sigma = np.std(Xelite ,axis=0,ddof=1).reshape(1,-1)
    
    # smoothing the mean
    xmean = Smoothing(xmean, xmean0, alpha)
    
    # dynamic smoothing parameter
    beta_t = beta * (1 - (1 - 1/t)**q)
    
    # smoothing the standard deviation
    sigma = Smoothing(sigma, sigma0, beta_t)
    
    return xmean, xmedian, xbest, Fmean, Fmedian, Fbest, sigma
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# Smoothing - smoothing scheme for variable update
# -----------------------------------------------------------------
def Smoothing(xnew, xold, s):

    # apply a smoothing scheme based on the parameter s
    return s * xnew + (1 - s) * xold
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# ComputeErrorS - compute standard deviation error
# -----------------------------------------------------------------
def ComputeErrorS(sigma, sigma0, TolAbs, TolRel):

    # error weights vector
    ewt = ErrorWeights(sigma,TolAbs,TolRel)
    
    # standard deviation error
    ErrorS = wrmsNorm(sigma - sigma0, ewt)
    
    # convergence metric based on standard deviation
    SmallErrorS = ErrorS <= 1.0
    
    return ErrorS, SmallErrorS
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# ErrorWeights - vector with the error weights
# -----------------------------------------------------------------
def ErrorWeights(x, TolAbs, TolRel):

    return 1.0/(TolAbs+np.abs(x)*TolRel)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# wrmsNorm - weighted root-mean-square norm
# -----------------------------------------------------------------
def wrmsNorm(v, w):

    return np.linalg.norm(v * w) / np.sqrt(v.size)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# UpdateLagrangeMult - update Lagrange multipliers
# -----------------------------------------------------------------
def UpdateLagrangeMult(G, H, lambdaG, lambdaH, Penalty):
        
    # update Lagrange multipliers for equality constraints
    lambdaH = lambdaH + Penalty*H
    
    # update Lagrange multipliers for inequality constraints 
    lambdaG = np.maximum(0, lambdaG + Penalty*G)
    
    return lambdaG, lambdaH
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# ComputeErrorC - compute constraint error 
# -----------------------------------------------------------------
def ComputeErrorC(G, H, lambdaG, lambdaH, Penalty, TolCon, ErrorC0):
    
    # constraints violation metrics
    ViolationEqNorm = np.max(np.abs(H))
    ViolationInNorm = np.max(np.abs(np.minimum(-G,lambdaG/Penalty)))
    #ViolationEqNorm = np.linalg.norm(H, np.inf)
    #ViolationInNorm = np.linalg.norm(np.minimum(-G, lambdaG / Penalty), np.inf)
    
    # constraints violation error
    ErrorC = max(ViolationEqNorm,ViolationInNorm)
    
    # convergence indicator for constraint violation
    SmallErrorC = ErrorC <= TolCon*ErrorC0
    
    return ErrorC, SmallErrorC
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# UpdatePenalty - update penalty parameter
# -----------------------------------------------------------------
def UpdatePenalty(Penalty, PenaltyFactor, MaximumPenalty, SmallErrorC):

    if not SmallErrorC:
        Penalty = min(PenaltyFactor*Penalty,MaximumPenalty)
        
    return Penalty
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# CheckConv - verify the convergence
# -----------------------------------------------------------------
def CheckConv(Fopt, SmallErrorS, SmallErrorC, CEstr):

    if SmallErrorC is None:
        SmallErrorC = True
    
    ExitFlag = 0
    
    if CEstr['iter'] >= CEstr['MaxIter']:
        ExitFlag = 1
        return ExitFlag
    if CEstr['stall'] >= CEstr['MaxStall']:
        ExitFlag = 2
        return ExitFlag
    if CEstr['Fcount'] >= CEstr['MaxFcount']:
        ExitFlag = 3
        return ExitFlag
    if CEstr['iter'] >= CEstr['MaxStall']:
        Idx1 = CEstr['iter']
        Idx0 = CEstr['iter'] - CEstr['MaxStall'] + 1
        if np.ptp(CEstr['Fbest'][Idx0:Idx1+1, 0]) <= CEstr['TolFun'] and SmallErrorC:
            ExitFlag = 4
            return ExitFlag
    if SmallErrorS and SmallErrorC:
        ExitFlag = 5
        return ExitFlag
    if Fopt <= CEstr['MinFval']:
        ExitFlag = 6
        return ExitFlag
    
    return ExitFlag
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# PrintProgress - print iteration progress on the screen
# -----------------------------------------------------------------
def PrintProgress(t, Nvars, CEstr):

    # print header in the first level
    if t == 1 and Nvars <= 5:
        if CEstr['isConstrained']:
            print("\n iter   best f(x)       std dev variat   error constr    design variable(s) \n")
        else:
            print("\n iter   best f(x)       std dev variat   design variable(s) \n")
    elif t == 1 and Nvars > 5:
        print("It is not possible to print more than 5 design variables on the screen")
        if CEstr['isConstrained']:
            print("\n iter   best f(x)       std dev variat   error constr")
        else:
            print("\n iter   best f(x)       std dev variat")
    
    # initial string with (t, F, Error)
    if CEstr['isConstrained']:
        MyString = "\n {:5d} {:+.9E} {:.9E} {:.9E}"
    else:
        MyString = "\n {:5d} {:+.9E} {:.9E}"
    
    # print values on screen
    if Nvars <= 5:
        # append format for design variables
        for i in range(Nvars):
            MyString += " {:+.6E}"
        # values with x
        if CEstr['isConstrained']:
            print(MyString.format(t, CEstr['Fbest'][t-1, 0], CEstr['ErrorS'][t-1, 0], CEstr['ErrorC'][t-1, 0], *CEstr['xbest'][t-1, :]), end='')
        else:
            print(MyString.format(t, CEstr['Fbest'][t-1, 0], CEstr['ErrorS'][t-1, 0], *CEstr['xbest'][t-1, :]), end='')
    else:
        # values without x
        if CEstr['isConstrained']:
            print(MyString.format(t, CEstr['Fbest'][t-1, 0], CEstr['ErrorS'][t-1, 0], CEstr['ErrorC'][t-1, 0]), end='')
        else:
            print(MyString.format(t, CEstr['Fbest'][t-1, 0], CEstr['ErrorS'][t-1, 0]), end='')
# -----------------------------------------------------------------
            

# -----------------------------------------------------------------
# PrintEnd - Display a summary of the optimization results
# -----------------------------------------------------------------
def PrintEnd(Xopt, Fopt, ExitFlag, CEstr):

    if ExitFlag == 1:
        Msg = 'Maximum number of iterations reached. '
    elif ExitFlag == 2:
        Msg = 'Solution stalled: no significant change in objective function over a set number of iterations.'
    elif ExitFlag == 3:
        Msg = 'Maximum number of function evaluations reached.'
    elif ExitFlag == 4:
        Msg = 'Objective function range has not changed significantly after many iterations. '
        if CEstr['isConstrained']:
            Msg += 'Additionally, constraint violations are small, indicating a potential solution. '
    elif ExitFlag == 5:
        Msg = 'Standard deviation variation is small, suggesting convergence towards a solution. '
        if CEstr['isConstrained']:
            Msg += 'Constraint violations are also small, indicating satisfactory adherence to constraints. '
    elif ExitFlag == 6:
        Msg = 'Minimum function value criterion met. '
    else:
        Msg = 'Unknown termination reason. '
    print("\n\n" + Msg + "\n")
    
    if ExitFlag == 1 or ExitFlag == 3:
        print("\nConsider increasing the maximum number of iterations or function evaluations.")
    elif ExitFlag == 2 or ExitFlag == 4 or ExitFlag == 5:
        print("\nSolution appears to be optimal within specified tolerances.")
    elif ExitFlag == 6:
        print("\nOptimization successfully found a solution meeting the minimum function value criterion.")
    
    print("\n\n--------------------------------------------------------")
    print(" Summary of the Optimization Process with the CE method ")
    print("--------------------------------------------------------")
    
    # Display the optimal point found
    out_str = "Optimal Point  Found: {:+.6E} \n".format(Xopt[0,0])
    for i in range(Xopt.shape[1]-1):
        out_str += "                      {:+.6E} \n".format(Xopt[0,i+1])
    print(out_str)
    
    #out_str = "Optimal Point  Found: {:+.6E} \n".format(Xopt[0,0])
    #for i, xi in enumerate(np.ravel(Xopt)-1):
    #    out_str += "                      {:+.6E} \n".format(xi+1)
    #print(out_str)

    # Display the optimal value found
    print("Optimal Value  Found: {:+.6E}".format(float(Fopt)))
    print("\n")
 
    # Display the number of iterations performed
    print("Iterations Performed: {}".format(CEstr['iter'  ]))
    print("\n")
    
    # Display the number of stall iterations
    print("Iterations on  Stall: {}".format(CEstr['stall' ]))
    print("\n")
    
    # Display the total number of function evaluations
    print("Function Evaluations: {}".format(CEstr['Fcount']))
    print("--------------------------------------------------------")
# -----------------------------------------------------------------
    

# -----------------------------------------------------------------
# DeleteEmptyEntries - delete empty entries from sampling records
# -----------------------------------------------------------------
def DeleteEmptyEntries(t, CEstr):

    if t < CEstr['MaxIter']:
        CEstr['xmean'  ] = CEstr['xmean'  ][:t, :]
        CEstr['xmedian'] = CEstr['xmedian'][:t, :]
        CEstr['xbest'  ] = CEstr['xbest'  ][:t, :]
        CEstr['Fmean'  ] = CEstr['Fmean'  ][:t, :]
        CEstr['Fmedian'] = CEstr['Fmedian'][:t, :]
        CEstr['Fbest'  ] = CEstr['Fbest'  ][:t, :]
        CEstr['sigma'  ] = CEstr['sigma'  ][:t, :]
        CEstr['ErrorS' ] = CEstr['ErrorS' ][:t, :]
        CEstr['ErrorC'] = CEstr['ErrorC'][:t, :]
    if not CEstr['isConstrained']:
        CEstr['ErrorC'] = np.zeros((0, 1))
        
    return CEstr
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# DomainSampling - sample from truncated Gaussian distribution
# -----------------------------------------------------------------
def DomainSampling(mu, sigma, lb, ub, Nvars, Ns, X):

    # limit vectors for standard truncated Gaussian
    # Compute l and u as arrays of shape (Ns, Nvars)
    l = np.tile((lb - mu) / sigma, (Ns, 1))
    u = np.tile((ub - mu) / sigma, (Ns, 1))
    
    # generate samples from truncated Gaussian distribution
    for n in range(Nvars):
        X[:, n] = mu[0, n] + sigma[0, n]*trandn(l[:, n], u[:, n])
    
    return X
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# trandn
# -----------------------------------------------------------------
# This function is an efficient generator of a random vector of 
# dimension len(l) = len(u) from the standard multivariate
# normal distribution, truncated over the region [l,u]. Infinite
# values for bounds 'u' and 'l' are accepted.
# 
# Remark:
# If you wish to simulate a random variable Z from the 
# non-standard Gaussian N(m,s^2) conditional on l < Z < u, then 
# first simulate X = trandn((l-m)/s,(u-m)/s) and set Z = m + s*X.
# 
# Input:
# l - (Nvars x 1) lower bound (as 1D array)
# u - (Nvars x 1) upper bound (as 1D array)
# 
# Output:
# x - (Nvars x 1) random vector with multivariate distribution N(0,1)
# 
# References:
# Botev, Z. I. (2016). "The normal law under linear restrictions: 
# simulation and estimation via minimax tilting". Journal of the 
# Royal Statistical Society: Series B (Statistical Methodology). 
# https://doi.org/10.1111/rssb.12162
# 
# MATLAB Central File Exchange:
# Z. Botev, Truncated Normal Generator
# shorturl.at/hntuB
# -----------------------------------------------------------------
def trandn(l, u):

    #l = np.asarray(l).reshape(-1,1)
    #u = np.asarray(u).reshape(-1,1)
    if l.size != u.size:
        raise ValueError('Truncation limits have to be vectors of the same length')
    x = np.full(l.shape, np.nan)
    a = 0.66  # threshold for switching between methods
    # case 1: a < l < u
    I = l > a
    if np.any(I):
        tl   = l[I]
        tu   = u[I]
        x[I] = ntail(tl, tu)
    # case 2: l < u < -a
    J = u < -a
    if np.any(J):
        tl   = -u[J]
        tu   = -l[J]
        x[J] = -ntail(tl, tu)
    # case 3: otherwise use inverse transform or accept-reject
    I2 = ~(I | J)
    if np.any(I2):
        tl = l[I2]
        tu = u[I2]
        x[I2] = tn(tl, tu)
    return x

#  --- ntail --- 
# This function samples a column vector of dimension 
# len(l)=len(u) from the standard normal distribution, truncated over the region [l,u], where 
# l > 0 and l and u are arrays. It uses a sampling
# algorithm based on acceptance-rejection from a Rayleigh 
# distribution similar to Marsaglia (1964).
def ntail(l, u):
    c = l**2 / 2.0
    n = l.size
    f = np.expm1(c - (u**2)/2.0)
    x = c - np.log(1 + np.random.rand(n) * f)
    # keep list of rejected
    I = np.where((np.random.rand(n)**2 * x) > c)[0]
    while I.size > 0:
        cy        = c[I]
        y         = cy - np.log(1 + np.random.rand(I.size) * f[I])
        idx       = np.where((np.random.rand(I.size)**2 * y) < cy)[0]
        x[I[idx]] = y[idx]
        I         = np.delete(I, idx)

    x = np.sqrt(2 * x)

    return x

#  --- tn  --- 
# This function samples a column vector of dimension 
# len(l)=len(u) from the standard normal distribution, truncated over the region [l,u], where 
# -a < l < u < a for some 'a' and l and u are arrays.
# It uses acceptance-rejection and inverse-transform method.
def tn(l, u):

    tol = 2  # controls switch between methods
    x = np.copy(l)
    I = np.where(np.abs(u - l) > tol)[0]
    if I.size > 0:
        tl   = l[I]
        tu   = u[I]
        x[I] = trnd(tl, tu)
    I2 = np.where(np.abs(u - l) <= tol)[0]
    if I2.size > 0:
        tl    = l[I2]
        tu    = u[I2]
        pl    = erfc(tl / np.sqrt(2)) / 2.0
        pu    = erfc(tu / np.sqrt(2)) / 2.0
        x[I2] = np.sqrt(2) * erfcinv(2 * (pl - (pl - pu) * np.random.rand(I2.size)))

    return x

#  --- trnd ---
# This function uses an acceptance-rejection sampling strategy
# to simulate from a truncated normal.
def trnd(l, u):

    x = np.random.randn(*l.shape)
    I = np.where((x < l) | (x > u))[0]
    while I.size > 0:
        ly        = l[I]
        uy        = u[I]
        y         = np.random.randn(I.size)
        idx       = np.where((y > ly) & (y < uy))[0]
        x[I[idx]] = y[idx]
        I         = np.delete(I, idx)

    return x
# -----------------------------------------------------------------
