## Cross-Entropy Optimizer

**CEopt: Cross-Entropy Optimizer** is a Matlab package that implements a framework for nonconvex optimization using the Cross-Entropy (CE) method. Due to the algorithm's relative simplicity, CEopt provides a transparent "gray-box" optimization solver with intuitive control parameters. It effectively handles both equality and inequality constraints through an augmented Lagrangian method, offering robustness and scalability for moderately sized complex problems. **CEopt**'s applicability and effectiveness are demonstrated through select case studies, making it a practical addition to optimization research and application toolsets.

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

Customize the input parameters in the provided scripts to fit your specific optimization problem.

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
