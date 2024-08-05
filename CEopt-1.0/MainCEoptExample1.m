clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample1.m ')
disp(' ------------------- ')

% objective function
F = @(x) -0.8*exp(-(x-2).^2) - 0.5*exp(-(x+2).^2)+1;

% bound for design variables
lb = -5;
ub =  5;

% initialize mean and std. dev. vectors
mu0    = (ub+lb)/2;
sigma0 = (ub-lb)/6;

% CE optimizer
tic
[Xopt,Fopt,ExitFlag,CEobj] = CEopt(F,mu0,sigma0,lb,ub)
toc