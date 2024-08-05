% -----------------------------------------------------------------
%  MainCEoptExample2Ext.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Jun 18, 2021
%           Last updated in: Jul 24, 2024
% -----------------------------------------------------------------
%  ﻿Example 2: The Peaks function in 2D
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample2.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

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
CEstr.MaxIter      = 80;      % Maximum number of iterations
CEstr.TolAbs       = 1.0e-2;  % Absolute tolerance
CEstr.TolRel       = 1.0e-2;  % Relative tolerance
CEstr.alpha        = 0.7;     % Smoothing parameter
CEstr.beta         = 0.8;     % Smoothing parameter
CEstr.q            = 10;      % Smoothing parameter

% CE optimizer
tic
[Xopt, Fopt, ExitFlag, CEstr] = CEopt(F, mu0, sigma0, lb, ub, [], CEstr)
toc

% meshgrid for visualization
Ngrid  = 100;
xRange = linspace(lb(1),ub(1),Ngrid);
yRange = linspace(lb(2),ub(2),Ngrid);
[X, Y] = meshgrid(xRange, yRange);
Z = F([X(:), Y(:)]);
Z = reshape(Z, size(X));

% custom colors
MyYellow = [0.9290 0.6940 0.1250];
MyBlue   = [0.0000 0.4470 0.7410];
MyRed    = [0.6350 0.0780 0.1840];

% custom colormap
N = 50;
MyColorMap1 = viridis(N);

% set up the figure
figure;
hold on;

% plotting F(x)
[c1,fig1] = contourf(X, Y, Z, N);
colormap(MyColorMap1);
colorbar;

% highlighting the optimal point
fig2 = plot(0.22,-1.62,'o','Color',MyYellow,...
                          'MarkerSize', 15,...
                          'MarkerFaceColor',MyYellow,...
                          'DisplayName', 'Analytical Solution');

% plot initial approximation
fig3 = plot(mu0(1),mu0(2),'d', ...
                          'Color',MyRed,...
                          'LineWidth' , 1,...
                          'MarkerSize', 7,...
                          'MarkerFaceColor',MyRed,...
                          'DisplayName', 'Initial Approximation');

% plot CE iterations
for j=1:CEstr.iter
    fig4 = plot(CEstr.xmean(j,1), CEstr.xmean(j,2),'d', ...
                                                   'Color',MyRed, ...
                                                   'LineWidth',1,...
                                                   'MarkerSize',7, ...
                                                   'DisplayName','CE Iterations');
end

% plot CE solution
fig5 = plot(Xopt(1), Xopt(2),'x','Color',MyBlue,...
                             'LineWidth' , 3,...
                             'MarkerSize', 25,...
                             'DisplayName','CE Optimum');

% labeling
xlabel('x_1', 'FontSize', 20, 'FontName', 'Helvetica');
ylabel('x_2', 'FontSize', 20, 'FontName', 'Helvetica');

% legend
leg = [fig2; fig3; fig4; fig5];
leg = legend(leg,'Location','NorthEast','FontSize',12);

% Set font and box
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on;
hold off;

% save the figure
exportgraphics(gca, 'CEoptExample2.eps', 'Resolution', 300);

% Objective function
function F = PeaksFunc(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = 3*(1-x1).^2.*exp(-x1.^2 - (x2+1).^2) ...
      - 10*(x1/5 - x1.^3 - x2.^5).*exp(-x1.^2 - x2.^2)...
      - (1/3)*exp(-(x1+1).^2 - x2.^2);
end