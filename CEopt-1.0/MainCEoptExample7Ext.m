% -----------------------------------------------------------------
%  MainCEoptExample7Ext.m
% -----------------------------------------------------------------
%  programmer: Marcos Vinicius Issa
%              marcosviniciusissa@gmail.com
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Aug 23, 2024
% -----------------------------------------------------------------
%  ﻿Example 7: Nonconvex structural optimization
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample7.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% convertion factor
Inch2Meter = 0.0254;    % inch to meter factor

% truss parameters
L                  = 360.0*Inch2Meter; % bar length (m)
MyTruss.l1         = L;                % 1st length (m)
MyTruss.l2         = L;                % 2nd length (m)
MyTruss.h          = L;                % height (m)
MyTruss.rho        = 2770.0;           % material density (kg/m^3)
MyTruss.E          = 69.8e9;           % elastic modulus (Pa)
MyTruss.AddedMass  = 454.0;            % added mass (kg)
MyTruss.omegaTresh = 2*pi*[7 15 20]';  % treshold frequencies (rad/s)
MyTruss.FixedDoFs  = [9 10 11 12];     % Diriclet boundary condition DoFs
MyTruss.NODES      = [2*L,L;           % nodes coordinates
                      2*L,0; 
                      L,L; 
                      L,0; 
                      0,L; 
                      0,0];
MyTruss.ELEM       = [5,3;             % element nodes
                      3,1; 
                      6,4; 
                      4,2; 
                      4,3; 
                      2,1; 
                      5,4; 
                      3,6; 
                      3,2; 
                      1,4]; 
MyTruss.Nnodes      = size(MyTruss.NODES,1); % # of nodes
MyTruss.Nelem       = size(MyTruss.ELEM,1);  % # of elements
MyTruss.Ndofs       = 2*MyTruss.Nnodes;      % # of DoFs

% objective and constraint functions
fun     = @(x)TrussMass      (x,MyTruss);
nonlcon = @(x)TrussConstraint(x,MyTruss);

% number of variables
Nvars = 10;

% lower and upper bounds for design variables
lb    = 0.1*ones(1,Nvars)*Inch2Meter^2;
ub    =  70*ones(1,Nvars)*Inch2Meter^2;

% initial mean and standard deviation
xmean0 = 0.5*(ub+lb);
sigma0 = 5*(ub-lb);

% CE solver
tic
[Xopt,Fopt,ExitFlag,CEstr] = CEopt(fun,xmean0,sigma0,lb,ub,nonlcon);
toc

% check constraint violation
disp(' ')
disp('Check if inequality constraints are <= 0')
G = nonlcon(Xopt);
for i=1:length(G)
  fprintf('G(%d) =  %5.2f\n',i,G(i));
end

% plot the nominal truss structure
PlotTruss10(xmean0*400,MyTruss,'Non-optimal Truss Structure');

% plot the optimal truss structure
PlotTruss10(  Xopt*400,MyTruss,'Optimal Truss Structure');