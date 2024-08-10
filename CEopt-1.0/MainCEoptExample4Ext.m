% -----------------------------------------------------------------
%  MainCEoptExample4Ext.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Jun 18, 2021
%           Last updated in: Jul 24, 2024
% -----------------------------------------------------------------
%  ﻿Example 4: Identification of a harmonic oscillator
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample4.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% system response observations
wn = 1; ksi = 0.1; y0 = 2; v0 = 0.5;
t0 = 0.0; t1 = 30.0; Ndata = 50;
time = linspace(t0,t1,Ndata)';
w     = wn*sqrt(1-ksi^2);
A     = sqrt(y0^2 + ((v0+ksi*wn*y0)/w)^2);
phi   = atan(y0*w/(v0+ksi*wn*y0));
ytrue = A*exp(-ksi*wn*time).*sin(w*time+phi);
ydata = ytrue + 0.05*y0*randn(Ndata,1);

% objective function (minimize the discrepancy)
F = @(x) MyMisfitFunc(x,ydata,time);

% bound for design variables
lb = [0.0; 0.0; -3; -3];
ub = [2.0; 1.0;  3;  3];

% cross-entropy optimizer struct
CEstr.Nsamp  = 50;      % number of samples
CEstr.TolAbs = 1.0e-3;  % absolute tolerance
CEstr.TolRel = 1.0e-2;  % relative tolerance

% CE optimizer
tic
[Xopt,Fopt,ExitFlag,CEobj] = CEopt(F,[],[],lb,ub,[],CEstr)
toc

% animation of the parameters identification process
for n=1:CEobj.iter

    % model parameters
    wn     = CEobj.xmean(n,1);
    ksi    = CEobj.xmean(n,2);
    y0     = CEobj.xmean(n,3);
    v0     = CEobj.xmean(n,4);

    % ODE model
    dydt = @(t,y) [0 1; -wn^2 -2*ksi*wn]*y;

    % initial conditions
    IC = [y0 v0];

    % model response
    [time,ymodel] = ode45(dydt,time,IC);
    y             = ymodel(:,1);

    figure(1)
    plot(time,ytrue,'--k','LineWidth',2.0);
    hold on
    plot(time,ydata, 'or','LineWidth',2.0);
    plot(time,y    , '-b','LineWidth',2.0);
    hold off
    ylim([-3 3])
    leg = legend('real system','noisly data','identified model');
    set(leg,'FontSize',16);
    set(gca,'FontSize',18);
    saveFileName = sprintf('CEoptExample4_%d.eps',n);
    exportgraphics(gca, saveFileName, 'Resolution', 300);
end

% Misfit function
function J = MyMisfitFunc(x,ydata,tspan)
    [Ns,Nvars] = size(x);     % input dimensions
             J = zeros(Ns,1); % preallocate memory for J
            wn = x(:,1);      % model parameter 1
           ksi = x(:,2);      % model parameter 2
            y0 = x(:,3);      % model parameter 3
            v0 = x(:,4);      % model parameter 4
    for n = 1:Ns
        dydt = @(t,y) [0 1; -wn(n)^2 -2*ksi(n)*wn(n)]*y;
        [time,ymodel] = ode45(dydt,tspan,[y0(n) v0(n)]);
        J(n) = norm(ydata-ymodel(:,1))/sqrt(length(ydata));
    end
end