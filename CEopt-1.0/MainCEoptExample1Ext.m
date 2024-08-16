% -----------------------------------------------------------------
%  MainCEoptExample1Ext.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Jun 18, 2021
%           Last updated in: Jul 24, 2024
% -----------------------------------------------------------------
%  ﻿Example 1: A Gaussian mixture in 1D
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample1.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

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
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(F,mu0,sigma0,lb,ub)
toc

% domain for plotting
xgrid = linspace(lb,ub,1000);
ygrid = F(xgrid);

% color map
N = 15;
cmap = magma(N);

% set up the figure
figure
hold on

% plotting the Gaussian distributions
for i = 1:4:N
    mu      = CEstr.xmean(i);
    sigma   = CEstr.sigma(i);
    y_gauss = exp(-(xgrid-mu).^2/(2*sigma^2)) / sqrt(2*pi*sigma^2);
    plot(xgrid, y_gauss, 'Color', cmap(N+1-i,:), 'LineWidth', 0.8);
end

% plotting F(x)
plot(xgrid, ygrid, 'b-', 'LineWidth', 2);

% highlighting the optimal point
plot(Xopt, Fopt,'ko','MarkerFaceColor','k','MarkerSize',8);

% labeling
xlabel('x'   ,'FontSize',20,'FontName','Helvetica');
ylabel('F(x)','FontSize',20,'FontName','Helvetica');

% setting plot limits
xlim([lb ub]);
%ylim([0 max(ygrid)+1]);
ylim([0 25]);

% set font and box
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on
hold off

% save figure
exportgraphics(gca, 'CEoptExample1.eps', 'Resolution', 300);