clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample5.m ')
disp(' ------------------- ')

% data poitns
load Example5data.mat
xE     = X(:,1)';
yE     = X(:,2)';
thetaA = 0:360/(length(xE)-1):360;

% objective function
F = @(x)MyObjFunc(x,[xE; yE; thetaA]);

% bound for design variables
lb = [1.3; 2.0; 120.0; -15.0; -50.0];
ub = [2.5; 3.5; 150.0; - 5.0; -35.0];

% CE optimizer
tic
[Xopt, Fopt, ExitFlag, CEstr] = CEopt(F, [], [], lb, ub);
toc
