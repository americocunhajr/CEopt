% -----------------------------------------------------------------
%  MainCEoptExample8Ext.m
% -----------------------------------------------------------------
%  programmer: Julio Cesar de Castro Basilio
%              juliobasilio@ymail.com
%
%  Originally programmed in: Jul 13, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  ﻿﻿Example 8: Fractional-order controller optimal tuning
%  (this example may take several minutes to run)
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample8.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% system parameters
global Kp Ki Kd OrderI OrderD
m       = 0.1;
M       = 2.0;
L       = 0.5;
J       = 0.006;
r       = 0.0;
xdisp0  = 0.0;
dxdisp0 = 0.0;
theta0  = 10.0*pi/180;
dtheta0 =  0.0*pi/180;

% simulation window
tfinal  = 5.0;

% auxiliar parameters
A = (M+m)*(J+m*L^2);
B = m*L;
C = M+m;

% objective and constraint functions
fun     = @(x) MyObjFunc(x);
nonlcon = @(x) MyConstFunc(x);

% lower and upper bounds for design variables
lb = [-100; -100; -100; 0.001; 0.001];
ub = [ 100;  100;  100; 1.990; 1.990];

% cross-entropy optimizer struct
CEstr.EliteFactor  = 0.1;     % Elite samples percentage
CEstr.Nsamp        = 50;      % Number of samples
CEstr.MaxIter      = 80;      % Maximum number of iterations
CEstr.TolAbs       = 1.0e-2;  % Absolute tolerance
CEstr.TolRel       = 1.0e-2;  % Relative tolerance

% CE solver
tic
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(fun,[],[],lb,ub,nonlcon,CEstr);
toc

% Plotting the results
figure;
subplot(2,1,1);
plot(0:, theta, 'b-', 'LineWidth', 1.5);
hold on;
plot(tout, r, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angle \theta (rad)');
legend('\theta', '\theta_{ref}');
title('Angle \theta vs. Time');
grid on;

subplot(2,1,2);
plot(tout, xdisp, 'b-', 'LineWidth', 1.5);
hold on;
plot(tout, x_ref, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Displacement x (m)');
legend('x', 'x_{ref}');
title('Displacement x vs. Time');
grid on;

% objective function
function F = MyObjFunc(x)
    global Kp Ki Kd OrderI OrderD
    Kp     = x(1);
    Ki     = x(2);
    Kd     = x(3);
    OrderI = x(4);
    OrderD = x(5);
    
    sim('FOPIDControllerCartPendulum.slx');
     
    F = trapz(tout,e.^2) + trapz(tout,(usignal/100).^2);
end

% inequality constraints
function [G,H] = MyConstFunc(x)
    global Kp Ki Kd OrderI OrderD
    Kp     = x(1);
    Ki     = x(2);
    Kd     = x(3);
    OrderI = x(4);
    OrderD = x(5);
    
    sim('FOPIDControllerCartPendulum.slx');
    
    G = max(abs(xdisp)) - 0.5;
    H = [];
end