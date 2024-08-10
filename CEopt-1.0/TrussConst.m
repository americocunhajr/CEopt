% -----------------------------------------------------------------
%  Truss10Const.m
% -----------------------------------------------------------------
%  programmer: Marcos Vinicius Issa
%              marcos.issa@uerj.br
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  Constraints for truss optimization problem.
% -----------------------------------------------------------------
function [G,H] = TrussConst(A,MyTruss)

    % truss structure parameters
    rho       = MyTruss.rho;
    E         = MyTruss.E;
    AddedMass = MyTruss.AddedMass;
    NODES     = MyTruss.NODES;
    ELEM      = MyTruss.ELEM; 
    Nnodes    = MyTruss.Nnodes;
    Nelem     = MyTruss.Nelem;
    Ndofs     = MyTruss.Ndofs;

    % initialize inequatily constraint
    G = zeros(1,3);
    
    % initialize equatily constraint
    H = zeros(1,Ndofs);
    
    % preallocate memory for global matrices
    K = zeros(Ndofs,Ndofs); 
    M = zeros(Ndofs,Ndofs);
    
    % assembly global matrices
    for e = 1:Nelem
        % distance between NODESs
        dx = NODES(ELEM(e,2),1) - NODES(ELEM(e,1),1);
        dy = NODES(ELEM(e,2),2) - NODES(ELEM(e,1),2);
        l  = sqrt(dx^2+dy^2);
        % strain matrix
        c = dx/l; 
        s = dy/l;
        B = 1/l*[-c -s c s];
        % element DoFs
        eDof = [2*ELEM(e,1)-1, 2*ELEM(e,1),...
                2*ELEM(e,2)-1, 2*ELEM(e,2)];
        % local stiffness and mass matrices
        Ke = B'*E*B*A(e)*l;
        Me = rho*l*A(e)*[2 0 1 0;
                         0 2 0 1;
                         1 0 2 0;
                         0 1 0 2]/6;
        % assembly global stiffness and mass matrices
        K(eDof,eDof) = K(eDof,eDof) + Ke;
        M(eDof,eDof) = M(eDof,eDof) + Me;
    end
    
    % update mass matrix with the added mass
    M = M + AddedMass*eye(Ndofs,Ndofs);
    
    % fixed and free nodes
    FIX  = union(2*[5 6]-1,2*[5 6]);
    FREE = setdiff(1:2*Nnodes,FIX);
    
    % solve the generalized eigenvalue problem
    [V,D] = eig(K(FREE,FREE),M(FREE,FREE));
    
    % sort frequencies
    [omega2,Idx] = sort(diag(D));

    % equilibrium constraints
    H = 0;
    
    % frequency constraints
    G(1) = 1 - sqrt(omega2(1))/(2*pi*7);
    G(2) = 1 - sqrt(omega2(2))/(2*pi*15);
    G(3) = 1 - sqrt(omega2(3))/(2*pi*20);
% -----------------------------------------------------------------