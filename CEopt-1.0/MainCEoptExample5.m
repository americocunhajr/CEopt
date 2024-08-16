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
[Xopt, Fopt, ExitFlag, CEstr] = CEopt(F, [], [], lb, ub)
toc

% objective function
function F = MyObjFunc(x,p)

    xE     = p(1,:);
    yE     = p(2,:);
    thetaA = p(3,:);
    a0     = 1;
    b0     = x(1);
    c0     = x(1);
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
    
    F1 = sum((xE1 - real(xEc1)).^2 + (yE1 - real(yEc1)).^2);
    F2 = sum(abs(imag(xEc1).^2) + abs(imag(yEc1).^2) + abs(imag(thetaB).^2));
    F  = F1 + F2;
end