% -----------------------------------------------------------------
%  CEopt.m
% -----------------------------------------------------------------
%  Programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Jun 18, 2021
%           Last updated in: Dec 11, 2024
% -----------------------------------------------------------------
%  This routine employs the Cross-entropy (CE) method to solve the 
%  following optimization problem:
%  
%                        Xopt = arg min F(x)
%   
%                        subject to
%  
%                        lb <= x <= ub
%                        g(x) <= 0
%                        h(x)  = 0
%  where:
%  - F   : R^Nvars -> R is a given scalar objective function
%  - lb  : (1 x Nvars) vector of lower bounds for the decision variables
%  - ub  : (1 x Nvars) vector of upper bounds for the decision variables
%  - g(x): R^Nvars -> R^m is a vector of m inequality constraint functions
%  - h(x): R^Nvars -> R^p is a vector of p equality constraint functions
%
%  The  goal  is  to  minimize the scalar objective function F(x) defined 
%  within  a  known  rectangular  domain  (feasible  region),  while also 
%  satisfying the given equality and inequality constraints. The feasible 
%  region  is  further  defined by the constraints g(x) <= 0 and h(x) = 0, 
%  in addition to the bounds lb <= x <= ub.
%
%  The  algorithm  samples the feasible region using a truncated Gaussian
%  distribution  and  updates its parameters with the aid of an elite set
%  (defined  by  the  better  samples), aiming to transform this Gaussian
%  distribution into  a Dirac distribution centered at the global optimum.
%  For  constrained  optimization  problems,  the  CE method integrates a 
%  mechanism  for  handling  the  constraints,  such   as  the  augmented 
%  Lagrangian  method, to incorporate the effects of constraints into the 
%  optimization process.
%  
%  Input:
%  fun    -  Function handle for  the  objective function. This function
%            must accept a 1 x Nvars row  vector  (representing a single 
%            sample) or an M x Nvars matrix (representing M samples with 
%            variables in columns) as input and return a scalar value or 
%            a row vector of M scalar values (for vectorized operations) 
%            respectively.
%  xmean0  - (1 x Nvars) initial mean
%  sigma0  - (1 x Nvars) intial standard deviation
%  lb      - (1 x Nvars) lower bound
%  ub      - (1 x Nvars) upper bound
%  nonlcon - nonlinear constraint function
%  CEstr   -  Struct containing parameters and settings for the CE method.
%  
%  CEstr fields include:
%  * Verbose          : boolean flag to enable/disable screem output
%  * isConstrained    : boolean flag to indicate a constrained problem
%  * isVectorized     : boolean flag to indicate a vectorized function
%  * Nvars            : number of design variables in x
%  * EliteFactor      : proportion of samples for the elite set
%  * Nsamp            : number of samples to draw per iteration
%  * MaxIter          : maximum number of iterations
%  * MaxStall         : maximum number of stall iterations
%  * MaxFcount        : maximum number of objective function evaluations
%  * MinFval          : minimum admissible value for objective function 
%  * TolAbs           : absolute tolerance
%  * TolRel           : relative tolerance
%  * TolCon           : constraint violation tolerance
%  * TolFun           : function value change tolerance
%  * alpha            : smoothing parameter for the mean update
%  * beta             : smoothing parameter for the std. dev. update
%  * q                : exponent in dynamic smoothing parameter
%  * NonlconAlgorithm : algorithm to handle nonlinear constraints
%  * InitialPenalty   : initial penalty value for the augmented Lagrangian
%  * PenaltyFactor    : factor by which the penalty parameter is increased
%  * MaximumPenalty   : maximum value for the penalty parameter
%  * xmean            : history of mean value over iterations
%  * xmedian          : history of median over iterations
%  * xbest            : history of best sample point over iterations
%  * Fmean            : history of objective function mean over elite set
%  * Fmedian          : history of objective function median over elite set
%  * Fbest            : history of objective function best value found
%  * sigma            : history of standard deviation over iterations
%  * ErrorS           : history of standard deviation error
%  * ErrorC           : history of constraint violation error (for constrained problems)
%  * iter             : total number of iterations performed
%  * stall            : number of iterations without significant progress
%  * Fcount           : total number of function evaluations performed
%  * ConvergenceStatus: boolean flag indicating if the algorithm converged
%                       (true) or not (false).
%            
%  This struct can also include additional fields for customized behavior 
%  or extensions of the CE method.
%  
%  Output:
%  Xopt     - (1 x Nvars) optimal point
%  Fopt     - scalar optimal value
%  ExitFlag - Flag indicating the reason for algorithm termination:
%             0 - algorithm is running or has not been initialized
%             1 - maximum number of iterations reached
%             2 - no significant change in objective function (stalled)
%             3 - maximum number of function evaluations reached
%             4 - function change and constranint error smaller than tolerance
%             5 - std dev and constraint errors smaller than tolerance
%             6 - minimum function value criterion met
%  CEstr   -  The updated Cross-Entropy object struct containing the final
%             state  of  the  algorithm and possibly additional diagnostic
%             information.
% -----------------------------------------------------------------
%  References:
%  
%  [1] Reuven Y. Rubinstein, Dirk P. Kroese,
%      The Cross-Entropy Method: A Unified Approach to Combinatorial 
%      Optimization, Monte-Carlo Simulation, and Machine Learning,
%      Springer-Verlag, 2004.
%  
%  [2] A. Cunha Jr, M. V. Issa, J. C. Basilio, J. G. Telles Ribeiro,
%      CEopt: A MATLAB Package for Non-convex Optimization with the
%      Cross-Entropy Method, ArXiv, 2024
% -----------------------------------------------------------------
%  Copyright (C) 2024  Americo Cunha Jr et al.
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program. If not, see <https://www.gnu.org/licenses/>.
% -----------------------------------------------------------------
function [Xopt,Fopt,ExitFlag,CEstr] = ...
                    CEopt(fun,xmean0,sigma0,lb,ub,nonlcon,CEstr)

    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 7
        error('Too many inputs.')
    end

    % consistency check for the mandatory inputs parameters
    [lb,ub,xmean0,sigma0,Nvars] = CheckInput(lb,ub,xmean0,sigma0);
    
    % check if CEstr is not provided or is empty
    if nargin < 7 || isempty(CEstr)
        CEstr = struct();
    end

    % set the default parameters for CEstr (if necessary)
    CEstr = InitializeCEstr(CEstr,Nvars);

    % consistency check for CEstr parameters
    CheckCEstr(CEstr);

    % check if nonlcon is not provided or is empty
    if nargin < 6 || isempty(nonlcon)
        CEstr.isConstrained = false;

        % objective function
        ObjFun = @(x) fun(x);
    else
        CEstr.isConstrained = true;

        % objective function
        if strcmp(CEstr.NonlconAlgorithm,'AugLagrLog')
            ObjFun = @(x,Lg,Lh,p) AugLagrangian1(x,fun,nonlcon,Lg,Lh,p);
        else 
            ObjFun = @(x,Lg,Lh,p) AugLagrangian2(x,fun,nonlcon,Lg,Lh,p);
        end
    end

    % decide the appropriate optimization solver
    if ~CEstr.isConstrained
        
        % CE solver for an unconstrained problem
        [Xopt,Fopt,ExitFlag,CEstr] = ...
         UnconstrSolverCE(ObjFun,Nvars,xmean0,sigma0,lb,ub,CEstr);
    else
        
        % CE solver for a constrained problem
        [Xopt,Fopt,ExitFlag,CEstr] = ...
         ConstrSolverCE(ObjFun,Nvars,xmean0,sigma0,lb,ub,nonlcon,CEstr);
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% CheckInput - verify input parameters for possibles errors
% -----------------------------------------------------------------
function [lb,ub,xmean0,sigma0,Nvars] = ...
                           CheckInput(lb,ub,xmean0,sigma0)

    % ensure lb and ub are row vectors
    lb = lb(:)';
    ub = ub(:)';
    
    % check if lb and ub are empty
    if isempty(lb) || isempty(ub)
        error('lb and ub must be non-empty')
    end
    
    % number of variables
    Nvars = length(lb);
    
    % check for consistency in lb and ub
    if Nvars ~= length(ub)
        error('lb and ub must have the same dimension')
    end
    if any(isnan(lb)) || any(isnan(ub))
        error('lb and ub cannot have a NaN components')
    end
    if any(lb >= ub)
        error('lb < ub for all components')
    end
    
    % define xmean0 (if necessary)
    if isempty(xmean0)
        xmean0                = (ub+lb)/2;
        xmean0(isinf(xmean0)) = 0.0;
        xmean0(isnan(xmean0)) = 0.0;
    end
    
    % define sigma0 (if necessary)
    if isempty(sigma0)
        sigma0                = (ub-lb)/sqrt(12);
        sigma0(isinf(sigma0)) = 1.0;
        sigma0(isnan(sigma0)) = 1.0;
    end
    
    % ensure xmean0 and sigma0 are row vectors
    xmean0 = xmean0(:)';
    sigma0 = sigma0(:)';

    % check for consistency in xmean0 and sigma0
    if size(xmean0) ~= size(sigma0)
        error('xmean0 and sigma0 must have the same dimensions')
    end
    if any(isnan(xmean0)) || any(isnan(sigma0))
        error('xmean0 and sigma0 cannot have a NaN components')
    end
    if any(isinf(xmean0)) || any(isinf(sigma0))
        error('xmean0 and sigma0 cannot have a Inf components')
    end
    if any(xmean0 < lb) || any(xmean0 > ub)
        error('xmean0 must be in [lb,ub] interval')
    end
    if any(sigma0 <= 0.0)
        error('All components of sigma0 must be positive')
    end
    
    % check for dimension consistency in xmean0
    if length(xmean0) ~= Nvars
        error('xmean0 must be a 1 x Nvars vector')
    end
    
    % check for dimension consistency in sigma0
    if length(sigma0) ~= Nvars
        error('sigma0 must be a 1 x Nvars vector')
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% InitializeCEstr - initialize and set default parameters for CEstr
% -----------------------------------------------------------------
function CEstr = InitializeCEstr(CEstr,Nvars)

    % names and values for default fields of CEstr
    DefaultParams = struct(     ...
        'Verbose'         ,true       , ...
        'isConstrained'   ,false      , ...
        'isVectorized'    ,false      , ...
        'Nvars'           ,Nvars      , ...
        'EliteFactor'     ,0.05       , ...
        'Nsamp'           ,100        , ...
        'MaxIter'         ,100*Nvars  , ...
        'MaxStall'        ,50         , ...
        'MaxFcount'       ,+Inf       , ...
        'MinFval'         ,-Inf       , ...
        'TolAbs'          ,1.0e-6     , ...
        'TolRel'          ,1.0e-3     , ...
        'TolCon'          ,1.0e-3     , ...
        'TolFun'          ,1.0e-3     , ...
        'alpha'           ,0.4        , ...
        'beta'            ,0.4        , ...
        'q'               ,10         , ...
        'NonlconAlgorithm','AugLagLog', ... 
        'InitialPenalty'  ,10         , ...
        'PenaltyFactor'   ,100        , ...
        'MaximumPenalty'  ,+Inf         ...
    );
    
    % assign values for undefined fields in CEstr
    DefaultFields = fieldnames(DefaultParams);
    for i = 1:numel(DefaultFields)
        if ~isfield(CEstr, DefaultFields{i})
            CEstr.(DefaultFields{i}) = DefaultParams.(DefaultFields{i});
        end
    end

    % get the strange fields in CEstr
    StrangeFields = rmfield(CEstr,DefaultFields);

    % remove the strange fields from CEstr
    CEstr = rmfield(CEstr,fieldnames(StrangeFields));

    % order the default fields in CEstr
    CEstr = orderfields(CEstr,DefaultParams);

    % preallocate memory or initiate objects
    CEstr.xmean       = NaN*ones(CEstr.MaxIter,Nvars);
    CEstr.xmedian     = NaN*ones(CEstr.MaxIter,Nvars);
    CEstr.xbest       = NaN*ones(CEstr.MaxIter,Nvars);
    CEstr.Fmean       = NaN*ones(CEstr.MaxIter,1    );
    CEstr.Fmedian     = NaN*ones(CEstr.MaxIter,1    );
    CEstr.Fbest       = NaN*ones(CEstr.MaxIter,1    );
    CEstr.sigma       = NaN*ones(CEstr.MaxIter,Nvars);
    CEstr.ErrorS      = NaN*ones(CEstr.MaxIter,1    );
    CEstr.ErrorC      = NaN*ones(CEstr.MaxIter,1    );
    CEstr.iter        = 0;
    CEstr.stall       = 0;
    CEstr.Fcount      = 0;
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% CheckCEstr - check parameters consistency for CEstr
% -----------------------------------------------------------------
function CheckCEstr(CEstr)

    % flag to show output on screem
    if mod(CEstr.Verbose,1) ~= 0
        error('Verbose must be integer')
    end

    % elite samples percentage
    if ~isnumeric(CEstr.PenaltyFactor)   || ...
                  CEstr.EliteFactor <= 0 || ...
                  CEstr.EliteFactor > 1
        error('EliteFactor must be such that 0 < EliteFactor <= 1')
    end

    % number of samples
    if ~isnumeric(CEstr.Nsamp) || CEstr.Nsamp <= 1
        error('Nsamp must be greather than 1')
    end

    % maximum number of iterations
    if mod(CEstr.MaxIter,1) ~= 0 || CEstr.MaxIter < 1
        error('MaxIter must be a positive integer')
    end

    % maximum number of stall iterations
    if mod(CEstr.MaxStall,1) ~= 0 || CEstr.MaxStall < 1
        error('MaxStall must be a positive integer')
    end

    % maximum number of function evaluations
    if (mod(CEstr.MaxFcount,1) == 0 && CEstr.MaxFcount < 1) || ...
       (mod(CEstr.MaxFcount,1) ~= 0 && CEstr.MaxFcount ~= Inf)
        error('MaxFcount must be a positive integer or infinity')
    end

    % minimum function value
    if ~isnumeric(CEstr.MinFval)
        error('MinFval must be numeric')
    end

    % absolute tolerance
    if ~isnumeric(CEstr.TolAbs) || any(CEstr.TolAbs <= 0.0)
        error('TolAbs must be positive real')
    end
    
    % relative tolerance
    if ~isnumeric(CEstr.TolRel) || CEstr.TolRel < 0.0
        error('TolRel must be non-negative real')
    end

    % constraint tolerance
    if ~isnumeric(CEstr.TolCon) || CEstr.TolCon < 0.0
        error('TolCon must be non-negative real')
    end

    % function tolerance
    if ~isnumeric(CEstr.TolFun) || CEstr.TolFun < 0.0
        error('TolFun must be non-negative')
    end

    % smoothing parameter (0 < alpha <= 1)
    if ~isnumeric(CEstr.alpha) || CEstr.alpha <= 0 || CEstr.alpha > 1
        error('alpha must be such that 0 < alpha <= 1')
    end

    % dynamic smoothing parameter
    if ~isnumeric(CEstr.beta) || CEstr.beta <= 0
        error('beta must be non-negative')
    end

    % dynamic smoothing parameter
    if ~isnumeric(CEstr.q) || CEstr.q <= 0
        error('q must be non-negative')
    end

    % nonlinear constraint algorithm 
    if  ~ischar(CEstr.NonlconAlgorithm)             || ...
       (~strcmp(CEstr.NonlconAlgorithm,'AugLagLog') && ...
        ~strcmp(CEstr.NonlconAlgorithm,'AugLagMax'))
        error('Unknown option for NonlconAlgorithm')
    end

    % initial penalty
    if ~isnumeric(CEstr.InitialPenalty) || CEstr.InitialPenalty <= 0
        error('InitialPenalty must be must be non-negative')
    end

    % penalty factor
    if ~isnumeric(CEstr.PenaltyFactor) || CEstr.PenaltyFactor <= 1
        error('PenaltyFactor must be greather than 1')
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% UnconstrSolverCE - solve an unconstrained optimization problem
% -----------------------------------------------------------------
function [Xopt,Fopt,ExitFlag,CEstr] = ...
              UnconstrSolverCE(fun,Nvars,xmean0,sigma0,lb,ub,CEstr)

    % initialize parameters
    t           = 0;                        % iteration counter
    stall       = 0;                        % stall iterations counter
    Fcount      = 0;                        % function evaluation counter
    EliteFactor = CEstr.EliteFactor;        % elite factor
    Nsamp       = CEstr.Nsamp;              % total number of samples
    Nelite      = round(EliteFactor*Nsamp); % number of elite samples
    MaxIter     = CEstr.MaxIter;            % maximum of iterations
    TolAbs      = CEstr.TolAbs;             % absolute   tolerance
    TolRel      = CEstr.TolRel;             % relative   tolerance
    alpha       = CEstr.alpha;              % smoothing parameter for mean
    beta        = CEstr.beta;               % smoothing parameter for std. dev.
    q           = CEstr.q;                  % dynamic update parameter
    Xopt        = NaN*xmean0;               % optimal point
    Fopt        = +Inf;                     % optimal value
    ExitFlag    = 0;                        % termination condition flag
    
    % preallocate memory for design variables samples
    X = zeros(Nsamp,Nvars);
    
    % preallocate memory for objective function (if necessary)
    if ~CEstr.isVectorized
        F = NaN*ones(Nsamp,1);
    end

    % loop to sample the domain and update the distribution
    while ExitFlag == 0 && t <= MaxIter

        % update level counter
        t = t + 1;

        % sample the domain from a truncated Gaussian distribution
        X = DomainSampling(xmean0,sigma0,lb,ub,Nvars,Nsamp,X);
        
        % evaluate objective function
        if ~CEstr.isVectorized
            % case where fun is not a vectorized function
            for n = 1:Nsamp
                F(n) = fun(X(n,:));
            end
        else
            % case where fun is a vectorized function
            F = fun(X);
        end
        
        % update function evalution counter
        Fcount = Fcount + Nsamp;
        
        % define elite samples set 
        EliteSetId = DefineEliteSet(F,Nelite);

        % update the distribution parameters
        [xmean,xmedian,xbest,Fmean,Fmedian,Fbest,sigma] = ...
         UpdateDistribution(F,X,EliteSetId,xmean0,sigma0,alpha,beta,q,t);
        
        % update standard deviation error
        [ErrorS,SmallErrorS] = ComputeErrorS(sigma,sigma0,TolAbs,TolRel);
        
        % update old parameters
        xmean0 = xmean;
        sigma0 = sigma;
        
        % update the optimum
        if Fbest < Fopt
           CEstr.xbest(t,:) = xbest;
           CEstr.Fbest(t)   = Fbest;
           Xopt             = xbest;
           Fopt             = Fbest;
           stall            = 0;
        else
            CEstr.xbest(t,:) = CEstr.xbest(t-1,:);
            CEstr.Fbest(t)   = CEstr.Fbest(t-1);
            stall            = stall + 1;
        end
        
        % update optimization process history
        CEstr.iter         = t;
        CEstr.stall        = stall;
        CEstr.Fcount       = Fcount;
        CEstr.xmean(t,:)   = xmean;
        CEstr.xmedian(t,:) = xmedian;
        CEstr.Fmean(t)     = Fmean;
        CEstr.Fmedian(t)   = Fmedian;
        CEstr.sigma(t,:)   = sigma;
        CEstr.ErrorS(t)    = ErrorS;

        % print iteration progress on the screen
        if CEstr.Verbose
            PrintProgress(t,Nvars,CEstr);
        end
        
        % check the convergence
        ExitFlag = CheckConv(Fopt,SmallErrorS,[],CEstr);
    end

    % convergence check and update of 'ConvergenceStatus' field
    if ExitFlag > 3
        CEstr.ConvergenceStatus = true;
    else
        CEstr.ConvergenceStatus = false;
    end
    
    % print resume
	if CEstr.Verbose
        PrintEnd(Xopt,Fopt,ExitFlag,CEstr);
	end
    
	% delete empty entries from sampling records
    CEstr = DeleteEmptyEntries(t,CEstr);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% ConstrSolverCE - solve a constrained optimization problem
% -----------------------------------------------------------------
function [Xopt,Fopt,ExitFlag,CEstr] = ...
           ConstrSolverCE(fun,Nvars,xmean0,sigma0,lb,ub,nonlcon,CEstr)

    % initialize parameters
    t           = 0;                 % iteration counter
    stall       = 0;                 % stall iterations counter
    Fcount      = 0;                 % function evaluation counter
    EliteFactor = CEstr.EliteFactor; % elite factor
    Nsamp       = CEstr.Nsamp;       % total number of samples
    Nelite      = EliteFactor*Nsamp; % number of elite samples
    MaxIter     = CEstr.MaxIter;     % maximum of iterations
    TolAbs      = CEstr.TolAbs;      % absolute   tolerance
    TolRel      = CEstr.TolRel;      % relative   tolerance
    TolCon      = CEstr.TolCon;      % constraint tolerance
    alpha       = CEstr.alpha;       % smoothing parameter for mean
    beta        = CEstr.beta;        % smoothing parameter for std. dev.
    q           = CEstr.q;           % dynamic update parameter
    Xopt        = NaN*xmean0;        % optimal point
    Fopt        = +Inf;              % optimal value
    ExitFlag    = 0;                 % termination condition flag
    
    % initialize penalty parameters
    Penalty        = CEstr.InitialPenalty;
    PenaltyFactor  = CEstr.PenaltyFactor;
    MaximumPenalty = CEstr.MaximumPenalty;

    % initialize Lagrange multipliers
    [G0,H0] = nonlcon(xmean0);
    if isempty(G0), G0 = 0.0; end
    if isempty(H0), H0 = 0.0; end
    lambdaG = zeros(size(G0));
    lambdaH = zeros(size(H0));

    % initialize constraint error
    ErrorC = ComputeErrorC(G0,H0,lambdaG,lambdaH,Penalty,TolCon,1.0);
    
    % preallocate memory for design variables samples
    X = zeros(Nsamp,Nvars);

    % preallocate memory for objective function (if necessary)
    if ~CEstr.isVectorized
        F = NaN*ones(Nsamp,1);
    end

    % loop to sample the domain and update the distribution
    while ExitFlag == 0 && t <= MaxIter

        % update level counter
        t = t + 1;

        % sample the domain from a truncated Gaussian distribution
        X = DomainSampling(xmean0,sigma0,lb,ub,Nvars,Nsamp,X);

        % evaluate augmented Lagrangian 
        if ~CEstr.isVectorized
            % case where fun is not a vectorized function
            for n = 1:Nsamp
                F(n) = fun(X(n,:),lambdaG,lambdaH,Penalty);
            end
        else
            % case where fun is a vectorized function
            F = fun(X,lambdaG,lambdaH,Penalty);
        end

        % update function evalution counter
        Fcount = Fcount + Nsamp;
        
        % define elite samples set 
        EliteSetId = DefineEliteSet(F,Nelite);

        % update the distribution parameters
        [xmean,xmedian,xbest,Fmean,Fmedian,Fbest,sigma] = ...
         UpdateDistribution(F,X,EliteSetId,xmean0,sigma0,alpha,beta,q,t);
        
        % standard deviation error
        [ErrorS,SmallErrorS] = ComputeErrorS(sigma,sigma0,TolAbs,TolRel);

        % evalute the constraints at xbest
        [G,H] = nonlcon(xbest);

        % evaluate the constraint error
        [ErrorC,SmallErrorC] = ...
                ComputeErrorC(G,H,lambdaG,lambdaH,Penalty,TolCon,ErrorC);
        
        % update Lagrange multipliers
        [lambdaG,lambdaH] = ...
                        UpdateLagrangeMult(G,H,lambdaG,lambdaH,Penalty);

        % update pentalty parameter
        Penalty = ...
          UpdatePenalty(Penalty,PenaltyFactor,MaximumPenalty,SmallErrorC);
        
        % update old parameters
        xmean0 = xmean;
        sigma0 = sigma;
        
        % update the optimum
        if Fbest < Fopt
           CEstr.xbest(t,:) = xbest;
           CEstr.Fbest(t)   = Fbest;
           Xopt             = xbest;
           Fopt             = Fbest;
           stall            = 0;
        else
            CEstr.xbest(t,:) = CEstr.xbest(t-1,:);
            CEstr.Fbest(t)   = CEstr.Fbest(t-1);
            stall            = stall + 1;
        end
        
        % update optimization process history
        CEstr.iter         = t;
        CEstr.stall        = stall;
        CEstr.Fcount       = Fcount;
        CEstr.xmean(t,:)   = xmean;
        CEstr.xmedian(t,:) = xmedian;
        CEstr.Fmean(t)     = Fmean;
        CEstr.Fmedian(t)   = Fmedian;
        CEstr.sigma(t,:)   = sigma;
        CEstr.ErrorS(t)    = ErrorS;
        CEstr.ErrorC(t)    = ErrorC;

        % print iteration progress on the screen
        if CEstr.Verbose
            PrintProgress(t,Nvars,CEstr);
        end
        
        % check the convergence
        ExitFlag = CheckConv(Fopt,SmallErrorS,SmallErrorC,CEstr);
    end

    % convergence check and update of 'ConvergenceStatus' field
    if ExitFlag > 3
        CEstr.ConvergenceStatus = true;
    else
        CEstr.ConvergenceStatus = false;
    end
    
    % print resume
	if CEstr.Verbose
        PrintEnd(Xopt,Fopt,ExitFlag,CEstr);
	end
    
	% delete empty entries from sampling records
    CEstr = DeleteEmptyEntries(t,CEstr);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% AugLagrangian1 - augmented Lagrangian function (log formulation)
% -----------------------------------------------------------------
function [AL,F,G,H] = ...
            AugLagrangian1(x,fun,nonlcon,lambdaG,lambdaH,Penalty)
    
    % objetive function
    F = fun(x);

    % nonlinear constraints
    [G,H] = nonlcon(x);
    if isempty(G), G = 0.0; end
    if isempty(H), H = 0.0; end

    % shift 
    s = lambdaG/Penalty;

    % augmented Lagrangian function
    AL = F - sum((s.*lambdaG).*log(s-G+eps),2) + ...
                             H*(lambdaH')+ 0.5*Penalty*sum(H.^2,2);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% AugLagrangian2 - augmented Lagrangian function (max formulation)
% -----------------------------------------------------------------
function [AL,F,G,H] = AugLagrangian2(x,fun,nonlcon,...
                                          lambdaG,lambdaH,Penalty)
    
    % objetive function
    F = fun(x);

    % nonlinear constraints
    [G,H] = nonlcon(x);
    if isempty(G), G = 0.0; end
    if isempty(H), H = 0.0; end

    % shifted equality constraint
    H_s = H + lambdaH/Penalty;

    % shifted inequality constraint
    G_s = G + lambdaG/Penalty;

    % augmented Lagrangian function
    AL = F + 0.5*Penalty*(sum(H_s.^2,2) + sum(max(0,G_s).^2,2));
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% DefineEliteSet - define elite samples set
% -----------------------------------------------------------------
function EliteSetId = DefineEliteSet(F,Nelite)

    % sort ObjFunc evaluations (order statistics)
    [Fsort,Isort] = sort(F);
    
    % elite samples indices
    EliteSetId = Isort(1:Nelite);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% UpdateDistribution - update distribution parameters
% -----------------------------------------------------------------
function [xmean,xmedian,xbest,Fmean,Fmedian,Fbest,sigma] = ...
          UpdateDistribution(F,X,EliteSetId,xmean0,sigma0,alpha,beta,q,t)
    
    % elite samples and values
    Felite = F(EliteSetId);
    Xelite = X(EliteSetId,:);

    % estimators for the ObjFunc minimum point
    xmean   =   mean(Xelite);
    xmedian = median(Xelite);
    xbest   = Xelite(1,:);

    % estimators for the ObjFunc minimum value
    Fmean   =   mean(Felite);
    Fmedian = median(Felite);
    Fbest   = Felite(1);

    % estimator for the standard deviation
    sigma = std(Xelite);

    % smoothing the mean
    xmean = Smoothing(xmean,xmean0,alpha);

    % dynamic smoothing parameter
    beta_t = beta*(1 - (1-1/t)^q);

    % smoothing the standard deviation
    sigma = Smoothing(sigma,sigma0,beta_t);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% Smoothing - smoothing scheme for variable update
% -----------------------------------------------------------------
function xnew = Smoothing(xnew,xold,s)
    
    % apply a smoothing scheme based on the parameter s
    xnew = s*xnew + (1-s)*xold;
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% ComputeErrorS - compute standard deviation error
% -----------------------------------------------------------------
function [ErrorS,SmallErrorS] = ...
                        ComputeErrorS(sigma,sigma0,TolAbs,TolRel)

    % error weigths vector
    ewt = ErrorWeights(sigma,TolAbs,TolRel);
        
    % standard devition error
    ErrorS = wrmsNorm(sigma-sigma0,ewt);

    % convergence metric based on standard deviation
    SmallErrorS = ErrorS <= 1.0;
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% ErrorWeights - vector with the error weights
% -----------------------------------------------------------------
function w = ErrorWeights(x,TolAbs,TolRel)

    % error weights vector
    w = 1./(TolAbs + abs(x)*TolRel);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% wrmsNorm - weigthed root-mean-square norm
% -----------------------------------------------------------------
function wnorm = wrmsNorm(v,w)
    
    % weighted-root-mean-square norm
    wnorm = norm(v.*w/length(v));
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% UpdateLagrangeMult - update Lagrange multipliers
% -----------------------------------------------------------------
function [lambdaG,lambdaH] = ...
                 UpdateLagrangeMult(G,H,lambdaG,lambdaH,Penalty)
    
    % check the nonlinear constraints
    if isempty(G), G = 0.0; end
    if isempty(H), H = 0.0; end

    % update Lagrange multipliers for equality constraints
    lambdaH = lambdaH + Penalty*H;
    
    % update Lagrange multipliers for inequality constraints 
    lambdaG = max(0,lambdaG + Penalty*G);
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% ComputeErrorC - compute constraint error
% -----------------------------------------------------------------
function [ErrorC,SmallErrorC] = ...
                 ComputeErrorC(G,H,lambdaG,lambdaH,Penalty,TolCon,ErrorC0)

    % check the nonlinear constraints
    if isempty(G), G = 0.0; end
    if isempty(H), H = 0.0; end

    % constraints violation metrics
    ViolationEqNorm = norm(H,Inf);
    ViolationInNorm = norm(min(-G,lambdaG/Penalty),Inf);
    
    % constraints violation error
    ErrorC = max(ViolationEqNorm,ViolationInNorm);

    % convergence indicator for constraint violation
    SmallErrorC = ErrorC <= TolCon*ErrorC0;
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% UpdatePenalty - update penalty parameter
% -----------------------------------------------------------------
function Penalty = ...
    UpdatePenalty(Penalty,PenaltyFactor,MaximumPenalty,SmallErrorC)

    % update penalty parameter
    if ~SmallErrorC
        Penalty = min(PenaltyFactor*Penalty,MaximumPenalty);
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% CheckConv - verify the convergence
% -----------------------------------------------------------------
function ExitFlag = CheckConv(Fopt,SmallErrorS,SmallErrorC,CEstr)

    % check if SmallErrorC is empty
    if isempty(SmallErrorC)
        SmallErrorC = true;
    end
    
    % default exit flag
    ExitFlag = 0;
    
    % check if maximum of iterations is reached
    if CEstr.iter >= CEstr.MaxIter
        ExitFlag = 1;
        return
    end

    % check if maximum of stall iterations is reached
    if CEstr.stall >= CEstr.MaxStall
        ExitFlag = 2;
        return
    end

    % check if the maximum of function evaluations is reached
    if CEstr.Fcount >= CEstr.MaxFcount
        ExitFlag = 3;
        return
    end

    % check if both the function change and constraint error are small
    if CEstr.iter >= CEstr.MaxStall
        
        Idx1 = CEstr.iter;
        Idx0 = CEstr.iter - CEstr.MaxStall + 1;

        if range(CEstr.Fbest(Idx0:Idx1)) <= CEstr.TolFun && SmallErrorC
            ExitFlag = 4;
            return
        end
    end

    % check if standard devition and constraint errors are small
    if SmallErrorS && SmallErrorC
        ExitFlag = 5;
        return
    end

    % check if the admissible minimum is reached
    if Fopt <= CEstr.MinFval
        ExitFlag = 6;
        return
    end
    
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% PrintProgress - print iteration progress on the screen
% -----------------------------------------------------------------
function PrintProgress(t,Nvars,CEstr)

    % print header in the first level
    if t == 1 && Nvars <= 5
        if CEstr.isConstrained
            MyHeader = ['\n iter   func value       error std dev',...
                        '   error constr    design variable(s) \n'];
        else
            MyHeader = ['\n iter   func value       error std dev',...
                        '   design variable(s) \n'];
        end
        fprintf(MyHeader);
    elseif t == 1 && Nvars > 5
        MyHeader = ['It is not possible to print more than',...
                            ' 5 design variables on the screen \n'];
        fprintf(MyHeader);
        if CEstr.isConstrained
            MyHeader = ['\n iter   func value       error std dev',...
                        '   error constr \n'];
        else
            MyHeader = '\n iter   func value       error std dev \n';
        end
        fprintf(MyHeader);
    end
    
    % initial string with (t, F, Error)
    if CEstr.isConstrained
        MyString = '\n %5g %+.9E %.9E %.9E';
    else
        MyString = '\n %5g %+.9E %.9E';
    end
    
	% print values on screen
    if Nvars <= 5
        % string with (t, F, err, x)
        for i=1:Nvars
            MyString = strcat(MyString,' %+.6E');
        end
        % values with x
        if CEstr.isConstrained
            fprintf(MyString,t,CEstr.Fmean(t) ,...
                               CEstr.ErrorS(t),...
                               CEstr.ErrorC(t),...
                               CEstr.xmean(t,:));
        else
            fprintf(MyString,t,CEstr.Fmean(t) ,...
                               CEstr.ErrorS(t),...
                               CEstr.xmean(t,:));
        end
    else
        % values without x
        if CEstr.isConstrained
            fprintf(MyString,t,CEstr.Fmean(t),CEstr.ErrorS(t),CEstr.ErrorC(t));
        else
            fprintf(MyString,t,CEstr.Fmean(t),CEstr.ErrorS(t));
        end
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% PrintEnd - Display a summary of the optimization results
% -----------------------------------------------------------------
function PrintEnd(Xopt, Fopt, ExitFlag, CEstr)
    
    % Interpret ExitFlag and display appropriate message
    switch ExitFlag
        case 1
            Msg = 'Maximum number of iterations reached. ';
        case 2
            Msg = ['Solution stalled: no significant change in ', ...
                   'objective function over a set number of iterations.'];
        case 3
            Msg = 'Maximum number of function evaluations reached.';
        case 4
            Msg = ['Objective function range has not changed significantly ', ...
                   'after many iterations. '];
            if CEstr.isConstrained
                Msg = [Msg, 'Additionally, constraint violations are small ,', ...
                            'indicating a potential solution. '];
            end
        case 5
            Msg = ['Standard deviation variation is small, suggesting ', ...
                   'convergence towards a solution. '];
            if CEstr.isConstrained
                Msg = [Msg,'Constraint violations are also small, ', ...
                           'indicating satisfactory adherence to constraints. '];
            end
        case 6
            Msg = 'Minimum function value criterion met. ';
        otherwise
            Msg = 'Unknown termination reason. ';
    end
    fprintf('\n\n%s\n',Msg);

    % Provide warnings or recommendations based on ExitFlag
    if ExitFlag == 1 || ExitFlag == 3
        Msg = ['\nConsider increasing the maximum number of iterations ',...
              'or function evaluations.'];
        fprintf(Msg);
    elseif ExitFlag == 2 || ExitFlag == 4 || ExitFlag == 5
        Msg = '\nSolution appears to be optimal within specified tolerances.';
        fprintf(Msg);
    elseif ExitFlag == 6
        Msg = ['\nOptimization successfully found a solution meeting ',...
               'the minimum function value criterion.'];
        fprintf(Msg);
    end

    disp(' ');
    disp(' ');
    disp('--------------------------------------------------------');
    disp(' Summary of the Optimization Process with the CE method ');
    disp('--------------------------------------------------------');
    
    % Display the optimal value found
    fprintf('Optimal Value  Found: %+.6E\n', Fopt);
    
    % Display the optimal point found
    fprintf('Optimal Point  Found: [');
    for i = 1:length(Xopt)
        fprintf('%+.6E ', Xopt(i));
        if mod(i, 5) == 0 && i ~= length(Xopt)
            % Break line for readability
            fprintf('...\n              ');
        end
    end
    fprintf(']\n');
    
    % Display the number of iterations performed
    fprintf('Iterations Performed: %d\n', CEstr.iter);

    % Display the number of stall iterations
    fprintf('Iterations on  Stall: %d\n', CEstr.stall);
    
    % Display the total number of function evaluations
    fprintf('Function Evaluations: %d\n', CEstr.Fcount);
    disp('--------------------------------------------------------');
end

% -----------------------------------------------------------------

% -----------------------------------------------------------------
% DeleteEmptyEntries - delete empty entries from sampling records
% -----------------------------------------------------------------
function CEstr = DeleteEmptyEntries(t,CEstr)
	
    if t < CEstr.MaxIter
        CEstr.xmean     = CEstr.xmean(1:t,:);
        CEstr.xmedian   = CEstr.xmedian(1:t,:);
        CEstr.xbest     = CEstr.xbest(1:t,:);
        CEstr.Fmean     = CEstr.Fmean(1:t,1);
        CEstr.Fmedian   = CEstr.Fmedian(1:t,1);
        CEstr.Fbest     = CEstr.Fbest(1:t,1);
        CEstr.sigma     = CEstr.sigma(1:t,:);
        CEstr.ErrorS    = CEstr.ErrorS(1:t,1);
        CEstr.ErrorC    = CEstr.ErrorC(1:t,1);
    end
    if ~CEstr.isConstrained
        CEstr.ErrorC = [];
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% DomainSampling - sample from truncated Gaussian distribution
% -----------------------------------------------------------------
function X = DomainSampling(mu,sigma,lb,ub,Nvars,Ns,X)
    
    % limit vectors for standard truncated Gaussian
    l = ones(Ns,1)*((lb - mu)./sigma);
    u = ones(Ns,1)*((ub - mu)./sigma);

    % generate samples from truncated Gaussian distribution
    for n=1:Nvars
        X(:,n) = mu(n) + sigma(n)*trandn(l(:,n),u(:,n));
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% trandn
% -----------------------------------------------------------------
%  This function is an efficient generator of a random vector of 
%  dimension length(l) = length(u) from the standard multivariate
%  normal distribution, truncated over the region [l,u]. Infinite
%  values for bounds 'u' and 'l' are accepted.
% 
%  Remark:
%  If  you  wish  to  simulate  a  random  variable  'Z'  from the 
%  non-standard Gaussian N(m,s^2) conditional  on  l < Z < u, then 
%  first simulate X = trandn((l-m)/s,(u-m)/s) and set Z = m + s*X.
% 
%  Input:
%  l - (Nvars x 1) lower bound
%  u - (Nvars x 1) upper bound
%  
%  Output:
%  x - (Nvars x 1) random vector with multiv. distribution N(0,1)
% 
%  References:
%  Botev, Z. I. (2016). "The normal law under linear restrictions: 
%  simulation and estimation via minimax tilting". Journal of the 
%  Royal Statistical Society: Series B (Statistical Methodology). 
%  https://doi.org/10.1111/rssb.12162
%  
%  MATLAB Central File Exchange:
%  Z. Botev, Truncated Normal Generator
%  shorturl.at/hntuB
% -----------------------------------------------------------------
function x = trandn(l,u)
    l = l(:); u = u(:); % make 'l' and 'u' column vectors
    if length(l)~=length(u)
        error('Truncation limits have to be vectors of the same length')
    end
    x = NaN(size(l));
    a = .66; % treshold for switching between methods
    % threshold can be tuned for maximum speed for each Matlab version
    % three cases to consider:
    % case 1: a < l < u
    I = l > a;
    if any(I)
        tl = l(I); tu = u(I); x(I) = ntail(tl,tu);
    end
    % case 2: l < u < -a
    J = u < -a;
    if any(J)
        tl=-u(J); tu = -l(J); x(J) = -ntail(tl,tu);
    end
    % case 3: otherwise use inverse transform or accept-reject
    I = ~(I|J);
    if any(I)
       tl = l(I); tu = u(I); x(I) = tn(tl,tu);
    end
end

%  --- ntail --- 
%  This function samples a column vector of dimension 
%  length=length(l)=length(u) from the standard multivariate
%  normal distribution, truncated over the region [l,u], where 
%  l > 0 and l and u are column vectors. It uses a sampling
%  algorithm based on acceptance-rejection from a Rayleigh 
%  distribution similar to Marsaglia (1964).
function x = ntail(l,u)
    c = l.^2/2; n = length(l); f = expm1(c-u.^2/2);
    x = c - reallog(1+rand(n,1).*f); % sample using Rayleigh
    % keep list of rejected
    I = find(rand(n,1).^2.*x>c); d = length(I);
    while d > 0           % while there are rejections
               cy = c(I); % find the thresholds of rejected
                y = cy - reallog(1+rand(d,1).*f(I));
              idx = rand(d,1).^2.*y<cy; % accepted
        x(I(idx)) = y(idx);             % store the accepted
                I = I(~idx);            % remove accepted from list
                d = length(I);          % number of rejected
    end
    x = sqrt(2*x); % this Rayleigh transform can be delayed till the end
end

% --- tn --- 
%  This function samples a column vector  of dimension 
%  length=length(l)=length(u) from the standard multivariate
%  normal distribution, truncated over the region [l,u], where 
%  -a < l < u < a for some 'a' and l and u are column vectors.
%  It uses acceptance rejection and inverse-transform method.
function x = tn(l,u)
    tol = 2; % controls switch between methods
    % threshold can be tuned for maximum speed for each platform
    % case: abs(u-l) > tol, uses accept-reject from randn
    I = abs(u-l) > tol; x = l;
    if any(I)
        tl = l(I); tu = u(I); x(I) = trnd(tl,tu);
    end
    % case: abs(u-l) < tol, uses inverse-transform
    I = ~I;
    if any(I)
          tl = l(I); tu = u(I); 
          pl = erfc(tl/sqrt(2))/2; pu = erfc(tu/sqrt(2))/2;
        x(I) = sqrt(2)*erfcinv(2*(pl-(pl-pu).*rand(size(tl))));
    end
end

% --- trnd ---
%  This function uses an acceptance-rejection sampling strategy
%  to simulate from truncated normal.
function x = trnd(l,u)
    x = randn(size(l)); % sample normal
    % keep list of rejected
    I = find(x < l | x > u); d=length(I);
    while d>0 % while there are rejections
               ly = l(I); % find the thresholds of rejected
               uy = u(I);
                y = randn(size(ly));
              idx = y > ly & y < uy; % accepted
        x(I(idx)) = y(idx);           % store the accepted
                I = I(~idx);          % remove accepted from list
                d = length(I);        % number of rejected
    end
end
% -----------------------------------------------------------------