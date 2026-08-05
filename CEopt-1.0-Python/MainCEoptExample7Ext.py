# -----------------------------------------------------------------
#  MainCEoptExample7Ext.py
# -----------------------------------------------------------------
#  programmer: Marcos Vinicius Issa
#              marcosviniciusissa@gmail.com
# 
#  Originally programmed in: Mar 28, 2025
#           Last updated in: Jun 12, 2025
# -----------------------------------------------------------------
#  Example 7: Nonconvex structural optimization
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.linalg import eig
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample7Ext.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
# Objective function
# -----------------------------------------------------------------
def TrussMass(A, MyTruss):
    
    A  = np.atleast_2d(A)
    
    rho   = MyTruss['rho'  ]
    NODES = MyTruss['NODES']
    ELEM  = MyTruss['ELEM' ]
    Nelem = MyTruss['Nelem']
    
    M_total = 0.0
    
    for e in range(Nelem):
        node1    = int(ELEM[e, 0]) - 1
        node2    = int(ELEM[e, 1]) - 1
        dx       = NODES[node2, 0] - NODES[node1, 0]
        dy       = NODES[node2, 1] - NODES[node1, 1]
        l        = np.sqrt(dx**2 + dy**2)
        M_total += rho * A[0,e] * l

    return M_total
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# Constraint function
# -----------------------------------------------------------------
def TrussConstraint(A, MyTruss):

    A = np.atleast_2d(A)
    A = np.ravel(A)

    rho        = MyTruss['rho'       ]
    E          = MyTruss['E'         ]
    AddedMass  = MyTruss['AddedMass' ]
    omegaTresh = MyTruss['omegaTresh']
    FixedDoFs  = MyTruss['FixedDoFs' ]
    NODES      = MyTruss['NODES'     ]
    ELEM       = MyTruss['ELEM'      ]
    Nelem      = MyTruss['Nelem'     ]
    Ndofs      = MyTruss['Ndofs'     ]
    
    Nconstr = len(omegaTresh)
    G = np.zeros((1,Nconstr))
    H = None  # no equality constraints
    
    K = np.zeros((Ndofs, Ndofs))
    M = np.zeros((Ndofs, Ndofs))
    
    for e in range(Nelem):
        node1 = int(ELEM[e, 0]) - 1
        node2 = int(ELEM[e, 1]) - 1
        dx    = NODES[node2, 0] - NODES[node1, 0]
        dy    = NODES[node2, 1] - NODES[node1, 1]
        l     = np.sqrt(dx**2 + dy**2)
        c     = dx / l
        s_val = dy / l
        B     = np.array([-c, -s_val, c, s_val])
        eDof  = [2*node1, 2*node1 + 1, 2*node2, 2*node2 + 1]
        Ke    = (E * A[e] / l) * np.outer(B, B)
        Me    = (rho * A[e] * l / 6) * np.array([[2, 0, 1, 0],
                                               [0, 2, 0, 1],
                                               [1, 0, 2, 0],
                                               [0, 1, 0, 2]])
        idx    = np.ix_(eDof, eDof)
        K[idx] += Ke
        M[idx] += Me

    M = M + AddedMass * np.eye(Ndofs)
    
    all_dofs = np.arange(1, Ndofs+1)
    FreeDoFs = np.setdiff1d(all_dofs, FixedDoFs) - 1  # convert to 0-indexed
    
    K_free  = K[np.ix_(FreeDoFs, FreeDoFs)]
    M_free  = M[np.ix_(FreeDoFs, FreeDoFs)]
    eigvals = eig(K_free, M_free, right=False)
    omega   = np.sort(np.sqrt(np.real(eigvals)))
    
    for j in range(Nconstr):
        G[0,j] = 1 - omega[j] / omegaTresh[j].item()
    
    return G, H
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# PlotTruss10 - Plot a 10-bar truss given the elements'
#               cross-section areas.
# -----------------------------------------------------------------
def PlotTruss10(x, MyTruss, MyTitle):
    
    x = np.asarray(x).ravel()

    # Unpack truss structure parameters
    h     = MyTruss['h'    ]
    NODES = MyTruss['NODES']
    ELEM  = MyTruss['ELEM' ]
    
    grayColor = [0.7, 0.7, 0.7]
    
    # Create a new figure with default axes font size 10
    fig, ax = plt.subplots(figsize=(8,6))
    ax.cla()
    
    # --- Plot Support 1 ---
    xl1 = [-1.0, -0.7]
    yl1 = [-1.0, 1.0]
    X1, Y1 = hatch_coordinates(xl1, yl1, 0.15)
    ax.plot(X1, Y1, 'k', linewidth=1.0)
    a11 = np.array([0, -0.7, -0.7])
    a12 = np.array([0,  1.0, -1.0])
    patch1 = patches.Polygon(np.column_stack((a11, a12)), closed=True,
                              facecolor=grayColor, edgecolor='none')
    ax.add_patch(patch1)
    ax.plot([a11[1], a11[2]], [a12[1], a12[2]], 'k', linewidth=2.5)
    
    # --- Plot Support 2 ---
    xl2 = [-1.0, -0.7]
    yl2 = [h - 1.0, h + 1.0]
    X2, Y2 = hatch_coordinates(xl2, yl2, 0.15)
    ax.plot(X2, Y2, 'k', linewidth=1.0)
    a21 = np.array([0, -0.7, -0.7])
    a22 = np.array([h, h + 1, h - 1])
    patch2 = patches.Polygon(np.column_stack((a21, a22)), closed=True,
                              facecolor=grayColor, edgecolor='none')
    ax.add_patch(patch2)
    ax.plot([a21[1], a21[2]], [a22[1], a22[2]], 'k', linewidth=2.5)
    
    # --- Plot each truss element ---
    # Assume that len(x) equals the number of elements.
    # If ELEM comes from Matlab (1-based indices), subtract 1.
    for i in range(len(x)):
        # Get the nodes for element i. Adjust indices if needed.
        elem_nodes = NODES[np.array(ELEM[i]) - 1, :]  # Use "-1" for 1-based indexing
        # Draw the element outline with line width proportional to x[i]
        patch_elem = patches.Polygon(elem_nodes, closed=True, fill=False,
                                     edgecolor=grayColor, linewidth=x[i])
        ax.add_patch(patch_elem)
        # Also plot markers at the element nodes
        ax.plot(elem_nodes[:, 0], elem_nodes[:, 1], 'o', markersize=10,
                markeredgecolor='blue', markerfacecolor='white', linestyle='None',
                markeredgewidth=3)
    
    # --- Plot supports at the first 4 nodes ---
    for i in range(4):
        ax.plot(NODES[i, 0], NODES[i, 1], 'o', markersize=17,
                markeredgecolor=[0.9290, 0.6940, 0.1250],
                markerfacecolor='none', linewidth=4.1)
    
    # Remove tick labels and spines
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.set_title(MyTitle, fontsize=18)
    plt.pause(1)  # Pause for 1 second

    return fig, ax
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# hatch_coordinates - Return coordinates for plotting a hatch 
#                     pattern. The pattern is created so that
#                     the ratio of xstep/ystep adjusts the line
#                     angle.
# -----------------------------------------------------------------
def hatch_coordinates(xlim, ylim, xstep=1, ystep=None, merge=True):

    if ystep is None:
        ystep = xstep
    xpos = np.arange(xlim[0], xlim[1] + xstep, xstep)
    ypos = np.arange(ylim[0], ylim[1] + ystep, ystep)
    nx = xpos.size
    ny = ypos.size
    # Create a base grid similar to the Matlab implementation.
    nanline = np.full(nx + ny - 3, np.nan)
    row1 = np.concatenate((np.full(ny-2, xpos[0]), xpos[:-1]))
    row2 = np.concatenate((ypos[::-1][1:ny-1], np.full(nx-1, ypos[-1])))
    # Create a second row (for a different edge) similarly:
    row1b = np.concatenate((xpos[1:], np.full(ny-2, xpos[-1])))
    row2b = np.concatenate((np.full(nx-1, ypos[0]), ypos[::-1][1:ny-1]))
    
    # Stack the rows and append the nanline as a separator.
    Xmat = np.vstack((row1, row1b, nanline))
    Ymat = np.vstack((row2, row2b, nanline))
    
    if merge:
        X = Xmat.flatten()
        Y = Ymat.flatten()
    else:
        X = Xmat
        Y = Ymat

    return X, Y
# -----------------------------------------------------------------

# Conversion factor
Inch2Meter = 0.0254  # inch to meter factor

# Truss parameters
L                     = 360.0 * Inch2Meter  # bar length in meters
MyTruss               = {}
MyTruss['l1'        ] = L
MyTruss['l2'        ] = L
MyTruss['h'         ] = L
MyTruss['rho'       ] = 2770.0  # material density (kg/m^3)
MyTruss['E'         ] = 69.8e9  # elastic modulus (Pa)
MyTruss['AddedMass' ] = 454.0   # added mass (kg)
# treshold frequencies (rad/s)
MyTruss['omegaTresh'] = 2 * np.pi * np.array([7, 15, 20]).reshape(-1, 1)
MyTruss['FixedDoFs' ] = np.array([9, 10, 11, 12])
# Node coordinates: 6 nodes with 2 coordinates each
MyTruss['NODES'     ] = np.array([
    [2*L,  L],
    [2*L,  0],
    [ L,   L],
    [ L,   0],
    [ 0,   L],
    [ 0,   0]
])
# Element connectivity: each row contains two node indices (Matlab indices, adjust if needed)
MyTruss['ELEM'       ] = np.array([
    [5, 3],
    [3, 1],
    [6, 4],
    [4, 2],
    [4, 3],
    [2, 1],
    [5, 4],
    [3, 6],
    [3, 2],
    [1, 4]
])
MyTruss['Nnodes'      ] = MyTruss['NODES'].shape[0]
MyTruss['Nelem'       ] = MyTruss['ELEM'].shape[0]
MyTruss['Ndofs'       ] = 2 * MyTruss['Nnodes']

# Objective and constraint functions
fun     = lambda x: TrussMass(x, MyTruss)
nonlcon = lambda x: TrussConstraint(x, MyTruss)

# Number of design variables
Nvars = 10

# Bounds for design variables
lb = 0.1 * np.ones(Nvars) * (Inch2Meter**2)
ub = 70  * np.ones(Nvars) * (Inch2Meter**2)

# Initial mean and standard deviation
xmean0 = 0.5 * (ub + lb)
sigma0 =   5 * (ub - lb)

# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(fun, xmean0, sigma0, lb, ub, nonlcon)
elapsedTime = time.time() - start_time

print("\n")
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsedTime))

# Check constraint violation
print("\n")
print("--------------------------------------------------------")
print("Check if inequality constraints are <= 0")
print("--------------------------------------------------------")
G, _ = nonlcon(Xopt)
for i, g_val in enumerate(np.ravel(G), start=1):
    print(f"G({i}) = {g_val:5.2f}")
print("--------------------------------------------------------")

# Plot the nominal (non-optimal) truss structure
fig1, ax1 = PlotTruss10(xmean0*400, MyTruss, "Non-optimal Truss Structure")

# Plot the optimal truss structure
fig2, ax2 = PlotTruss10(  Xopt*400, MyTruss, "Optimal Truss Structure")

# Display both figures at once
plt.show()
