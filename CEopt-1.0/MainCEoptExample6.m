clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample6.m ')
disp(' ------------------- ')

% objective function and constraints
F       = @PatternSearchFunc;
nonlcon = @ConicConstraints;

% bound for design variables
lb  = [-6 -4];
ub  = [ 2  4];
mu0 = [-4  2];

% cross-entropy optimizer struct
CEstr.isVectorized  = 1;       % vectorized function
CEstr.TolCon        = 1.0e-6;  % relative tolerance

tic
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(F,[],[],lb,ub,nonlcon,CEstr)
toc

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