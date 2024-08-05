% -----------------------------------------------------------------
%  MainCEoptExample5Ext.m
% -----------------------------------------------------------------
%  programmer: Jose Geraldo Telles Ribeiro
%              jose.gt.ribeiro@gmail.com
%
%  Originally programmed in: Jun 20, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  ﻿Example 5: Design of a complex mechanism
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample5.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproduc0b0lity)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% data points
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
[Xopt, Fopt, ExitFlag, CEstr] = CEopt(F, [], [], lb, ub)
toc

% plot optimal mechanism
a0     = 1;
b0     = Xopt(1)*a0;
c0     = b0;
d0     = Xopt(2)*a0;
gamma  = Xopt(3);
thetaD = Xopt(4);
theta0 = Xopt(5);

thetaA_1 = 0:0.5:360;
thetaE   = (180-gamma)/2;
e0       = 2*b0*cosd(thetaE);
m        = sqrt(a0^2+d0^2-2*a0*d0*cosd(thetaA_1+theta0));
beta     = asind(a0*sind(thetaA_1+theta0)./m);
thetaB   = acosd((b0^2+m.^2-c0^2)./(2*b0*m))-beta;

xEc0 = a0*cosd(thetaA_1+theta0+thetaD)+e0*cosd(thetaB+thetaE+thetaD);

a   = a0*(max(xE)-min(xE))/(max(xEc0)-min(xEc0));
b   = b0*(max(xE)-min(xE))/(max(xEc0)-min(xEc0));
c   = c0*(max(xE)-min(xE))/(max(xEc0)-min(xEc0));
d   = d0*(max(xE)-min(xE))/(max(xEc0)-min(xEc0));
e   = e0*(max(xE)-min(xE))/(max(xEc0)-min(xEc0));
xEc = a*cosd(thetaA_1+theta0+thetaD)+e*cosd(thetaB+thetaE+thetaD);
yEc = a*sind(thetaA_1+theta0+thetaD)+e*sind(thetaB+thetaE+thetaD);
x_A = mean(xE)-mean(xEc);
y_A = mean(yE)-mean(yEc);
xEc = xEc+x_A;
yEc = yEc+y_A;

% set up the figure
figure;
hold on;
plot(xE ,yE ,'bo','LineWidth',2.0);
plot(xEc,yEc,'r-','LineWidth',2.0);
xlim([0 50]); ylim([0 40]);

% labeling
xlabel('x', 'FontSize', 20, 'FontName', 'Helvetica');
ylabel('y', 'FontSize', 20, 'FontName', 'Helvetica');

% legend
leg = legend('Target Trajectory','Designed Trajectory');
set(leg,'FontSize',12);

% Set font and box
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 18);
box on;
hold off

% save the figure
exportgraphics(gca, 'CEoptExample5.eps', 'Resolution', 300);

% objective function
function F = MyObjFunc(x,p)

    xE     = p(1,:);
    yE     = p(2,:);
    thetaA = p(3,:);
    a0     = 1;
    b0     = x(1);
    c0     = b0;
    d0     = x(2);
    gamma  = x(3);
    thetaD = x(4);
    theta0 = x(5);
    
    thetaE = (180-gamma)/2;
    e0     = 2*b0*cosd(thetaE);
    m      = sqrt(a0^2 + d0^2 - 2*a0*d0*cosd(thetaA + theta0));
    beta   = asind(a0*sind(thetaA + theta0)./m);
    thetaB = acosd((b0^2 + m.^2 - c0^2)./(2*b0*m)) - beta;
    xEc0   = a0*cosd(thetaA+theta0+thetaD)+e0*cosd(thetaB+thetaE+thetaD);
    yEc0   = a0*sind(thetaA+theta0+thetaD)+e0*sind(thetaB+thetaE+thetaD);
    
    xE1  = (xE   - mean(xE)  )/(max(xE)   - min(xE)  );
    yE1  = (yE   - mean(yE)  )/(max(xE)   - min(xE)  );
    xEc1 = (xEc0 - mean(xEc0))/(max(xEc0) - min(xEc0));
    yEc1 = (yEc0 - mean(yEc0))/(max(xEc0) - min(xEc0));
    
    F = sum(abs(xE1 - xEc1) + abs(yE1 - yEc1));
end