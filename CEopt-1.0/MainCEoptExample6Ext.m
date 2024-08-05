% -----------------------------------------------------------------
%  MainCEoptExample6Ext.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Mar 21, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  ﻿Example 6: Nonsmooth function with conic constraints
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample6.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% bound for design variables and initial mean
lb  = [-6 -4];
ub  = [ 2  4];
mu0 = [-4  2];

% objective function and constraints
F       = @PatternSearchFunc;
nonlcon = @ConicConstraints;

% cross-entropy optimizer struct
CEstr.isVectorized  = 1;       % vectorized function
CEstr.TolCon        = 1.0e-6;  % constraint tolerance

tic
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(F,mu0,[],lb,ub,nonlcon,CEstr)
toc

% meshgrid for visualization
Ngrid  = 100;
xRange = linspace(lb(1),ub(1),Ngrid);
yRange = linspace(lb(2),ub(2),Ngrid);
[X,Y]  = meshgrid(xRange,yRange);
Z = F([X(:), Y(:)]);
Z = reshape(Z, size(X));
[C,Ceq] = nonlcon([X(:),Y(:)]);
C   = reshape(C  ,size(X));
Ceq = reshape(Ceq,size(X));

% custom colors
MyYellow = [0.9290 0.6940 0.1250];
MyBlue   = [0.0000 0.4470 0.7410];
MyRed    = [0.6350 0.0780 0.1840];

% custom colormap
N = 50;
MyColorMap1 = viridis(N);
MyColorMap2 = inferno(N);

% set up the figure
figure;
hold on;

% plotting F(x)
[c1,fig1] = contourf(X,Y,Z,N,'DisplayName','Objective Function');
colormap(MyColorMap1);
colorbar;

% plotting H(x)=0
[c2,fig2] = contour(X,Y,Ceq,[0 0],...
                     'Color','k',...
                     'LineWidth',2,...
                     'DisplayName','  Equality Constraint');

% plotting G(x)<=0
[c3,fig3] = contourf(X,Y,-C,[0 0],...
                     'Color',MyYellow,...
                     'FaceColor', [1 0.9 0.8],...
                     'LineWidth',2,...
                     'DisplayName','Inequality Constraint');


% plot initial approximation
fig4 = plot(mu0(1),mu0(2),'d', ...
                          'Color',MyRed,...
                          'LineWidth' , 1,...
                          'MarkerSize', 7,...
                          'MarkerFaceColor',MyRed,...
                          'DisplayName', 'Initial Approximation');

% plot CE iterations
for j=1:CEstr.iter
    fig5 = plot(CEstr.xmean(j,1), CEstr.xmean(j,2),'d', ...
                                                   'Color',MyRed, ...
                                                   'LineWidth',1,...
                                                   'MarkerSize',7, ...
                                                   'DisplayName','CE Iterations');
end

% plot CE solution
fig6 = plot3(Xopt(1),Xopt(2),Fopt,'x',...
                                  'Color',MyBlue,...
                                  'LineWidth' , 3,...
                                  'MarkerSize', 25,...
                                  'DisplayName','CE Optimum');

% labeling
xlabel('x_1', 'FontSize', 20, 'FontName', 'Helvetica');
ylabel('x_2', 'FontSize', 20, 'FontName', 'Helvetica');

% legend
leg = [fig2; fig3; fig4; fig5; fig6];
leg = legend(leg,'Location','SouthWest','FontSize',12);

% Set font and box
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on;
hold off;

% save the figure
exportgraphics(gca, 'CEoptExample6.eps', 'Resolution', 300);

% objective function
function F = PatternSearchFunc(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = zeros(size(x1,1),1);
    for i = 1:size(x,1)
        if  x1(i) < -5
            F(i) = (x1(i)+5).^2 + abs(x2(i));
        elseif x1(i) < -3
            F(i) = -2*sin(x1(i)) + abs(x2(i));
        elseif x1(i) < 0
            F(i) = 0.5*x1(i) + 2 + abs(x2(i));
        elseif x1 >= 0
            F(i) = .3*sqrt(x1(i)) + 5/2 + abs(x2(i));
        end
    end
end

% equality and inequality constraints
function [G,H] = ConicConstraints(x)
    x1 = x(:,1);
    x2 = x(:,2);
    G  = 2*x1.^2 + x2.^2 - 3;
    H  = (x1+1).^2 - (x2/2).^4;
end