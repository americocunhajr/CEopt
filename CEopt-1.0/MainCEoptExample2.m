clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample2.m ')
disp(' ------------------- ')

% objective function
F = @(x)PeaksFunc(x);

% bound for design variables
lb = [-3; -3];
ub = [ 3;  3];

% initialize mean and std. dev. vectors
mu0    = lb + (ub-lb).*rand(2,1);
sigma0 = 5*(ub-lb);
        
% define parameters for the CE optimizer
CEstr.isVectorized = 1;       % Vectorized function
CEstr.EliteFactor  = 0.1;     % Elite samples percentage
CEstr.Nsamp        = 50;      % Number of samples
CEstr.MaxIter      = 80;      % Maximum of iterations
CEstr.TolAbs       = 1.0e-2;  % Absolute tolerance
CEstr.TolRel       = 1.0e-2;  % Relative tolerance
CEstr.alpha        = 0.7;     % Smoothing parameter
CEstr.beta         = 0.8;     % Smoothing parameter
CEstr.q            = 10;      % Smoothing parameter

% CE optimizer
tic
[X_opt,F_opt,ExitFlag,CEstr] = CEopt(F,mu0,sigma0,lb,ub,[],CEstr)
toc

% objective function
function F = PeaksFunc(x)
    x1 = x(:,1);
    x2 = x(:,2);
     F = 3*(1-x1).^2.*exp(-x1.^2 - (x2+1).^2) ...
       - 10*(x1/5 - x1.^3 - x2.^5).*exp(-x1.^2 - x2.^2)...
       - (1/3)*exp(-(x1+1).^2 - x2.^2);
end