% -----------------------------------------------------------------
%  TrussMass.m
% -----------------------------------------------------------------
%  programmer: Marcos Vinicius Issa
%              marcos.issa@uerj.br
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Aug 23, 2024
% -----------------------------------------------------------------
%  Objective function for truss optimization problem.
% -----------------------------------------------------------------
function F = TrussMass(A,MyTruss)

    % truss structure parameters
    rho       = MyTruss.rho;
    NODES     = MyTruss.NODES;
    ELEM      = MyTruss.ELEM; 
    Nelem     = MyTruss.Nelem;
    
    % compute mass
    F = 0.0;
    for e = 1:Nelem
        dx = NODES(ELEM(e,2),1) - NODES(ELEM(e,1),1);
        dy = NODES(ELEM(e,2),2) - NODES(ELEM(e,1),2);
        l  = sqrt(dx^2+dy^2);
        F  = F + rho*A(e)*l;
    end
end
% -----------------------------------------------------------------