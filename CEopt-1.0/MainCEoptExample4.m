clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample4.m ')
disp(' ------------------- ')

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

% define parameters for the CE optimizer
CEstr.Nsamp  = 50;      % Number of samples
CEstr.TolAbs = 1.0e-3;  % Absolute tolerance
CEstr.TolRel = 1.0e-2;  % Relative tolerance

% CE optimizer
tic
[Xopt,Fopt,ExitFlag,CEobj] = CEopt(F,[],[],lb,ub,[],CEstr)
toc

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