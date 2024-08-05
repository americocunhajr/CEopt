% -----------------------------------------------------------------
%  MainCEoptExample9.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Jul 31, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  ﻿Example 9: SINDy
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample9.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% parameters and initial conditions
a  =  1.0;
b  = -1.0;
d  =  0.3;
f  =  0.65;
w  =  1;
IC = [1; 0; 0];

% temporal interval
t0 =   0.0;
t1 = 100.0;

% Forced Duffing oscillator
Duffing  = @(t,x)[x(2); -d*x(2)-b*x(1)-a*x(1)^3+f*cos(x(3)); w];

% integrate the ground truth
[t_true, X_true] = ode45(Duffing,[t0 t1],IC);

% generating training data
t_train = t_true(1:2:round((1/3)*end));
X_train = X_true(1:2:round((1/3)*end),:);
X       = X_train + 0.05*max(abs(X_train)).*randn(size(X_train));
dXdt    = zeros(size(X));
for n = 1:size(X,1)
    dXdt(n,:) = Duffing(0.0,X(n,:));
end
dXdt = dXdt + 0.05*max(abs(dXdt)).*randn(size(dXdt));

% sparsity promoting parameter
lambda = 0.25;

% dictionary of function
THETA_X = SINDyDictionary(X);

% objective function
J = @(z)MyMisfitFunc(z,THETA_X,dXdt,lambda);

% bounds for the parameters
XI_ls = reshape(THETA_X\dXdt,[],1);
lb    = XI_ls - 0.05*max(abs(XI_ls));
ub    = XI_ls + 0.05*max(abs(XI_ls));

% initial mean and standard deviation
mu0    = XI_ls;
sigma0 = (ub-lb)/sqrt(12);

% cross-entropy optimizer struct
CEstr.Verbose = 1;
CEstr.RelTol  = 1.0e-3;  % relative tolerance
CEstr.AbsTol  = 1.0e-6;  % relative tolerance
% CEstr.alpha   = 0.7;     % Smoothing parameter
% CEstr.beta    = 0.7;     % Smoothing parameter
% CEstr.q       = 10;      % Smoothing parameter

% CE solver
tic
XI_ce = CEopt(J, mu0, sigma0, lb, ub, [], CEstr);
toc

% number of base functions
Nbasis = size(THETA_X,2);

% number of state coordinates
Nstates = size(dXdt,2);

% identified coefficients via CE
fprintf('Identified coefficients:\n');
XI_ce = reshape(XI_ce,[Nbasis,Nstates])

% thresholded coefficients via CE
fprintf('Thresholded coefficients:\n');
small_coeff  = abs(XI_ce) < lambda;
XI_ce(small_coeff) = 0.0

% small coefficients in XI_ls
%       XI_ls = reshape(XI_ls,[],Nstates)
% small_coeff = abs(XI_ls) < lambda;
% XI_ls(small_coeff) = 0.0

% data-driven model
DataDrivenModel = @(t,x) (SINDyDictionary(x')*XI_ce)';

% integrate the identified model
[t_model,X_model] = ode45(DataDrivenModel,[t0 t1],IC);

% custom colors
MyYellow = [0.9290 0.6940 0.1250];
MyBlue   = [0.0000 0.4470 0.7410];
MyRed    = [0.6350 0.0780 0.1840];

% ploting x1(t)
figure(1)
hold on
plot(t_train,X_train(:,1),'o' ,'MarkerSize',7  ,'Color',MyYellow,'MarkerFaceColor',MyYellow);
plot(t_true ,X_true(:,1) ,'-' ,'LineWidth' ,1.2,'Color',MyBlue);
plot(t_model,X_model(:,1),'--','LineWidth' ,1.2,'Color',MyRed);
xlabel('time', 'FontSize', 20, 'FontName', 'Helvetica');
ylabel('x_1' , 'FontSize', 20, 'FontName', 'Helvetica');
legend({'training data','ground thuth','SINDy model'},'Location','NorthEast','FontSize',12);
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on
hold off
exportgraphics(gca, 'CEoptExample91.eps', 'Resolution', 300);

% ploting x2(t)
figure(2)
hold on
plot(t_train,X_train(:,2),'o' ,'MarkerSize',7  ,'Color',MyYellow,'MarkerFaceColor',MyYellow);
plot(t_true ,X_true(:,2) ,'-' ,'LineWidth' ,1.2,'Color',MyBlue);
plot(t_model,X_model(:,2),'--','LineWidth' ,1.2,'Color',MyRed);
xlabel('time', 'FontSize', 20, 'FontName', 'Helvetica');
ylabel('x_2' , 'FontSize', 20, 'FontName', 'Helvetica');
legend({'training data','ground thuth','SINDy model'},'Location','NorthEast','FontSize',12);
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on
hold off
exportgraphics(gca, 'CEoptExample92.eps', 'Resolution', 300);

% ploting x3(t)
figure(3)
hold on
plot(t_train,X_train(:,3),'o', 'MarkerSize',7  ,'Color',MyYellow,'MarkerFaceColor',MyYellow);
plot(t_true ,X_true(:,3) ,'-' ,'LineWidth' ,1.2,'Color',MyBlue);
plot(t_model,X_model(:,3),'--','LineWidth' ,1.2,'Color',MyRed);
xlabel('time', 'FontSize', 20, 'FontName', 'Helvetica');
ylabel('x_3' , 'FontSize', 20, 'FontName', 'Helvetica');
legend({'training data','ground thuth','SINDy model'},'Location','NorthEast','FontSize',12);
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on
hold off
exportgraphics(gca, 'CEoptExample93.eps', 'Resolution', 300);

% function to create the function dictionary
function THETA_X = SINDyDictionary(X)
    
    % number of state coordinates
    Nstates = size(X,2);

    % constant term
    THETA_X = [];
    THETA_X = [ones(size(X,1), 1), THETA_X];

    % polynomials up to the 3rd order
    for i = 1:Nstates
        THETA_X = [THETA_X, X(:,i)];
    end
    for i = 1:Nstates
        THETA_X = [THETA_X, X(:,i).^3];
    end
    
    % trigonometric terms
    THETA_X = [THETA_X, cos(X)];
end

% misfit function with zero-norm penalization
function J = MyMisfitFunc(z,THETA_X,dXdt,lambda)

    % number of state coordinates
    [Nt,Nstates] = size(dXdt);

    % matrix form of coefficients
    Z = reshape(z,[],Nstates);
    
    % mean-square error
    E2 = norm(dXdt-THETA_X*Z)/sqrt(Nt*Nstates);

    % zero-norm penalization
    ZeroNorm = sum(z > 0.0);
    
    % misfit function
    J = E2 + lambda*ZeroNorm;
end