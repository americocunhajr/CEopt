% -----------------------------------------------------------------
%  Truss10Mass.m
% -----------------------------------------------------------------
%  programmer: Marcos Vinicius Issa
%              marcos.issa@uerj.br
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  Objective function for truss optimization problem.
% -----------------------------------------------------------------
function M = TrussMass(A,MyTruss)

    % truss structure parameters
    rho       = MyTruss.rho;
    NODES     = MyTruss.NODES;
    ELEM      = MyTruss.ELEM; 
    Nelem     = MyTruss.Nelem;
    
    % compute mass
    M = 0.0;
    for e = 1:Nelem
        dx = NODES(ELEM(e,2),1) - NODES(ELEM(e,1),1);
        dy = NODES(ELEM(e,2),2) - NODES(ELEM(e,1),2);
        l  = sqrt(dx^2+dy^2);
        M = M + A(e)*l*rho;
    end
end
% -----------------------------------------------------------------