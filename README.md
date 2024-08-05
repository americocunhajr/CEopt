## Cross-Entropy Optimizer

**CEopt: Cross-Entropy Optimizer** is a Matlab package that implements a framework for nonconvex optimization using the Cross-Entropy (CE) method. Due to the algorithm's relative simplicity, CEopt provides a transparent "gray-box" optimization solver with intuitive control parameters. It effectively handles both equality and inequality constraints through an augmented Lagrangian method, offering robustness and scalability for moderately sized complex problems. **CEopt**'s applicability and effectiveness are demonstrated through select case studies, making it a practical addition to optimization research and application toolsets.

<p align="center">
<img src="logo/CEoptStructure.png" width="60%">
</p>

### Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Reproducibility](#reproducibility)
- [Documentation](#documentation)
- [Authors](#authors)
- [Citing CEopt](#citing-ceopt)
- [License](#license)
- [Institutional Support](#institutional-support)
- [Funding](#funding)
- [Contact](#contact)

### Overview
**CEopt** was developed to provide a robust and scalable solution for nonconvex optimization problems using the Cross-Entropy method. More details can be found in the following paper:
- **A. Cunha Jr, M. V. Issa, J. C. Basilio, and J. G. Telles Ribeiro**, *CEopt: A MATLAB Package for Nonconvex Optimization with the Cross-Entropy Method, 2024 (under review)*

Preprint available at: [Preprint Link](xxx)

### Features
- Implements Cross-Entropy method for nonconvex optimization
- Handles equality and inequality constraints using augmented Lagrangian method
- Transparent "gray-box" optimization solver with intuitive control parameters
- Robust and scalable for moderately sized complex problems
- Demonstrated applicability through select case studies

### Installation
To install and get started with **CEopt**, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/americocunhajr/CEopt.git
   ```
2. Navigate to the package directory:
   ```bash
   cd CEopt
   ```
3. Open Matlab and add the directory to the Matlab path:
   ```bash
   addpath(genpath('path_to_CEopt_directory/CEopt/CEopt-1.0'))
   ```

### Usage
To run CEopt, use the following commands in Matlab:
   ```bash
   [Xopt,Fopt,ExitFlag,CEstr] = CEopt(fun,xmean0,sigma0,lb,ub,nonlcon,CEstr)
   ```
Each parameter is described as follows:
- fun: Function handle for the objective function. This function must accept a 1 x Nvars row  vector  (representing a single sample) or an M x Nvars matrix (representing M samples with variables in columns) as input and return a scalar value or a row vector of M scalar values (for vectorized operations) respectively.
- xmean0: Initial mean of the design variables distributions.
- sigma0: Initial standard deviations for the design variables distributions.
- lb: Lower bounds for the design variables.
- ub: Upper bounds for the design variables.
- nonlcon: Function handle for the nonlinear constraints. Returns two arrays c (inequalities) and ceq (equalities).
- CEstr: Structure with settings for the CEopt algorithm.

The CEstr structure allows for extensive customization of the CE optimization process. Here's a breakdown of its fields:

- Verbose: Boolean flag to enable detailed output during optimization.
- isConstrained: Set to true if there are constraints defined by nonlcon.
- isVectorized: Set to true if fun and nonlcon can evaluate multiple rows of X in a single call.
- Nvars: Number of variables in the optimization problem.
- EliteFactor: Proportion of the population considered elite.
- Nsamp: Number of samples used in each iteration of the optimization.
- MaxIter: Maximum number of iterations for the optimization process.
- MaxStall: Maximum number of iterations without improvement before termination.
- MaxFcount: Maximum number of function evaluations allowed.
- MinFval: Target objective function value for early stopping.
- TolAbs: Absolute tolerance on the change in the objective function value for convergence.
- TolRel: Relative tolerance on the change in the objective function value for convergence.
- TolCon: Tolerance on the feasibility of constraints.
- TolFun: Tolerance on the change in function value for convergence.
- alpha: Smoothing parameter for the mean update.
- beta: Smoothing parameter for the standard deviation update.
- q: Smoothing parameter for standard deviation update.
- NonlconAlgorithm: Algorithm used for handling nonlinear constraints.
- InitialPenalty: Initial penalty coefficient for constraint violation.
- PenaltyFactor: Scaling factor for the penalty coefficient.

This extensive set of parameters and settings enables users to finely tune the CE optimization to their specific needs and problem characteristics.



### Reproducibility
The tutorials of **CEopt** package are fully reproducible. You can find a fully reproducible capsule of the simulations on <a href="https://codeocean.com/capsule/xxx" target="_blank">CodeOcean</a>.

### Documentation
The routines in **CEopt** package are well-commented to explain their functionality. Each routine includes a description of its purpose, as well as inputs and outputs. Detailed documentation about the code functionality can be found in paper.

### Authors
- Americo Cunha Jr
- Marcos Vinicius Issa
- Julio Cesar de Castro Basilio
- Jose Geraldo Telles Ribeiro

### Citing CEopt
If you use **CEopt** in your research, please cite the following publication:
- *A. Cunha Jr, M. V. Issa, J. C. Basilio and J. G. Telles Ribeiro, CEopt: A MATLAB Package for Nonconvex Optimization with the Cross-Entropy Method, 2024 (under review)*

```
@article{CunhaJr2024CEopt,
   author  = {A {Cunha~Jr} and M. V. Issa and J. C. Basilio and J. G. {Telles Ribeiro}},
   title   = "{CEopt: A MATLAB Package for Nonconvex Optimization with the Cross-Entropy Method}",
   journal = {Under Review},
   year    = {2024},
   volume  = {~},
   pages   = {~},
   doi    = {~},
}
```

### License
**CEopt** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

<img src="logo/mit_license_red.png" width="10%"> 

### Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%">

### Funding

<img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%">  &nbsp; &nbsp; &nbsp; <img src="logo/faperj.jpg" width="20%">

### Contact
For any questions or further information, please contact:

Americo Cunha Jr: americo.cunha@uerj.br
