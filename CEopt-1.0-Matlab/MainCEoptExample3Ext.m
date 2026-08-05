% -----------------------------------------------------------------
%  MainCEoptExample3Ext.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunhajr@gmail.com
%
%  Originally programmed in: Apr 01, 2024
%           Last updated in: Jul 24, 2024
% -----------------------------------------------------------------
%  ﻿Example 3: A collection of 2D algebraic benchmark problems
% -----------------------------------------------------------------

clc; clear; close all;

disp(' ------------------- ')
disp(' MainCEoptExample3.m ')
disp(' ------------------- ')

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
%rng_stream = RandStream('mt19937ar','Seed',30081985);
%rng_stream = RandStream('mt19937ar','Seed',30081986);
%rng_stream = RandStream('mt19937ar','Seed',30081987);
RandStream.setGlobalStream(rng_stream);

% List of benchmark functions, their bounds, and analytical optima
benchmarks = {
   {'Ackley'        , @Ackley        , [- 10.000,- 10.000], [+ 10.000,+ 10.000], {[+  0.00000,+  0.00000]},+  0.00000}
   {'Beale'         , @Beale         , [-  5.000,-  5.000], [+  5.000,+  5.000], {[+  3.00000,+  0.50000]},+  0.00000}
   {'Booth'         , @Booth         , [- 20.000,- 20.000], [+ 20.000,+ 20.000], {[+  1.00000,+  3.00000]},+  0.00000}
   {'BukinN6'       , @BukinN6       , [- 15.000,-  3.000], [-  5.000,+  3.000], {[- 10.00000,+  1.00000]},+  0.00000}
   {'CrossInTray'   , @CrossInTray   , [-  6.000,-  6.000], [+  6.000,+  6.000], {[-1.3491, -1.3491], [1.3491, 1.3491], [-1.3491, 1.3491], [1.3491, -1.3491]},-  2.06261}  % Multiple global optima
   {'DixonPrice'    , @DixonPrice    , [- 10.000,- 10.000], [+ 10.000,+ 10.000], {[+  1.00000,+  1.00000]},+  0.00000}
   {'Easom'         , @Easom         , [+  2.000,+  2.000], [+  4.000,+  4.000], {[+  3.14159,+  3.14159]},-  1.00000}
   {'Eggholder'     , @Eggholder     , [-512.000,-512.000], [+512.000,+512.000], {[+512.00000,+404.23190]},-959.64070}
   {'GoldsteinPrice', @GoldsteinPrice, [-  2.000,-  2.000], [+  2.000,+  2.000], {[+  0.00000,-  1.00000]},+  3.00000}
   {'Griewank'      , @Griewank      , [- 30.000,- 30.000], [+ 30.000,+ 30.000], {[+  0.00000,+  0.00000]},+  0.00000}
   {'Himmelblau'    , @Himmelblau    , [-  5.000,-  5.000], [+  5.000,+  5.000], {[+3.00000,+2.00000], [-2.805118,+3.283186], [-3.779310,-3.283186], [+3.584428, -1.848126]},+  0.00000}  % Multiple global optima
   {'HolderTable'   , @HolderTable   , [- 10.000,- 10.000], [+ 10.000,+ 10.000], {[+8.05502,+9.66459], [-8.055020,+9.664590], [+8.055020,-9.664590], [-8.055020, -9.664590]},- 19.20850}  % Multiple global optima
   {'LeviN13'       , @LeviN13       , [- 10.000,- 10.000], [+ 10.000,+ 10.000], {[+  1.00000,+  1.00000]},+  0.00000}
   {'Matyas'        , @Matyas        , [-100.000,-100.000], [+100.000,+100.000], {[+  0.00000,+  0.00000]},+  0.00000}
   {'McCormick'     , @McCormick     , [-  1.500,-  3.000], [+  4.000,+  4.000], {[-  0.54719,-  1.54719]},-  1.91330}
   {'Rastrigin'     , @Rastrigin     , [-  5.120,-  5.120], [+  5.120,+  5.120], {[+  0.00000,+  0.00000]},+  0.00000}
   {'Rosenbrock'    , @Rosenbrock    , [-  2.000,-  1.000], [+  2.000,+  3.000], {[+  1.00000,+  1.00000]},+  0.00000}
   {'SchafferN2'    , @SchafferN2    , [- 25.000,- 25.000], [+ 25.000,+ 25.000], {[+  0.00000,+  0.00000]},+  0.00000}
   {'SchafferN4'    , @SchafferN4    , [- 25.000,- 25.000], [+ 25.000,+ 25.000], {[+0.00000,+1.25313], [+0.00000,-1.25313], [+1.25313,+0.00000], [-1.25313,+0.00000]},+  0.29258}  % Approximate global optimum
   {'Shekel'        , @Shekel        , [-  2.000,-  2.000], [+ 12.000,+ 12.000], {[+  4.00000,+  4.00000]},- 10.15320}  % Approximate global optimum for m=10
   {'Sphere'        , @Sphere        , [-  5.000,-  5.000], [+  5.000,+  5.000], {[+  0.00000,+  0.00000]},+  0.00000}
   {'StyblinskiTang', @StyblinskiTang, [-  5.000,-  5.000], [+  5.000,+  5.000], {[-  2.90353,-  2.90353]},- 78.33230}  % For dimension n=2
   {'ThreeHumpCamel', @ThreeHumpCamel, [-  5.000,-  5.000], [+  5.000,+  5.000], {[+  0.00000,+  0.00000]},+  0.00000}
   {'Zakharov'      , @Zakharov      , [-  5.000,-  5.000], [+ 10.000,+ 10.000], {[+  0.00000,+  0.00000]},+  0.00000}
};

% Custom colormaps
MyColorMap1 = viridis(20);
%MyColorMap2 = inferno(20);

% Custom colors
MyYellow = [0.9290 0.6940 0.1250];
MyBlue   = [0.0000 0.4470 0.7410];
MyRed    = [0.6350 0.0780 0.1840];

% loop over the benchmark tests
for i = 1:length(benchmarks)

    % CEopt customization
    CEstr = struct();
    CEstr.Verbose = false;
    CEstr.Nsamp   = 1000;
    %CEstr.q = 20;
    %CEstr.alpha = 0.7;
    %CEstr.beta = 0.7;

    
    % Extract function info
    funcName           = benchmarks{i}{1};
    funcHandle         = benchmarks{i}{2};
    lb                 = benchmarks{i}{3};
    ub                 = benchmarks{i}{4};
    globalOptimaPoints = benchmarks{i}{5}; % global optima points
    globalOptimaValue  = benchmarks{i}{6}; % global optimum value
    
    % Number of variables 
    Nvars = length(lb);
    
    % Initialize mean and std. dev. vectors
    %mu0    = 0.5*(ub+lb);
    mu0    = lb + (ub-lb).*rand(1,Nvars);
    sigma0 = (ub-lb)/2;

    % Define an acceptable tolerance level for the optimization
    tol = 0.5; 
    
    % Display function being optimized
    disp(' '); 
    disp('---------------------------------------------------');
    disp([' Optimizing ', funcName, ' Function...']);
    
    % Call the CE optimizer
    tic
    [Xopt,Fopt,ExitFlag,CEstr] = ...
               CEopt(funcHandle,mu0,sigma0,lb,ub,[],CEstr);
    elapsedTime = toc;

    % Initialize minimum distance to infinity
    minDist = Inf;
    
    % Check if the solution is close to any known global optima
    closestOptimum = [];
    for opt = globalOptimaPoints
        dist = norm(Xopt - opt{1}); % distance to each global optimum
        if dist < minDist
            minDist = dist;
            closestOptimum = opt{1}; % Update the closest global optimum
        end
    end

    % Check if closest distance is within tolerance
    isCloseToOptimum = minDist < tol; 

    % Display results
    disp(' Numerical Solution:');
    disp([' Xopt: ', mat2str(Xopt)]);
    disp([' Fopt: ', num2str(Fopt)]);
    disp(' Analytical Solution:');
    for j = 1:length(globalOptimaPoints)
        disp([' X*: ', mat2str(globalOptimaPoints{j}),...
              ', F*: ', num2str(globalOptimaValue)]);
    end
    disp([' Closest global optimum: ', mat2str(closestOptimum)]);
    disp([' Minimum distance to global optimum: ', num2str(minDist)]);
    printExitFlagMeaning(ExitFlag);


    if isCloseToOptimum
        disp(' CEopt successfully approximated a global optimum.');
    else
        disp(' CEopt did not approximate a global optimum within the specified tolerance.');
    end
    disp('---------------------------------------------------');

    % Generate meshgrid for contour plot
    xRange = linspace(lb(1), ub(1), 100);
    yRange = linspace(lb(2), ub(2), 100);
    [X, Y] = meshgrid(xRange, yRange);
    Z = arrayfun(@(x,y) funcHandle([x y]), X, Y);
    
    figure; % Create new figure for each function
    hold on;
    
    % Contour plot of the function
    fig1 = contourf(X, Y, Z, 30); % Adjust the number of levels as needed
    colormap(MyColorMap1);
    colorbar;
    title([funcName, ' Function']);
    xlabel('x_1','FontSize',20,'FontName','Helvetica');
    ylabel('x_2','FontSize',20,'FontName','Helvetica');
    
    % Plot analytical global optima
    for opt = globalOptimaPoints
        fig2 = plot(opt{1}(1), opt{1}(2), 'o', ...
                               'Color',MyYellow,...
                               'MarkerSize', 15,...
                               'MarkerFaceColor',MyYellow,...
                               'DisplayName', 'Analytical Solution');
    end

    % Plot initial approximation
    fig3 = plot(mu0(1), mu0(2), 'd', ...
                                  'Color',MyRed,...
                                  'LineWidth' , 1,...
                                  'MarkerSize', 7,...
                                  'MarkerFaceColor',MyRed,...
                                  'DisplayName', 'Initial Approximation');

    % Plot CE iterations
    for j=1:CEstr.iter
        fig4 = plot(CEstr.xmean(j,1), CEstr.xmean(j,2), 'd', ...
                                  'Color',MyRed,...
                                  'LineWidth' , 1,...
                                  'MarkerSize', 7,...
                                  'DisplayName', 'CE Iterations');
    end

    % Plot CE solution
    fig5 = plot(Xopt(1), Xopt(2), 'x', ...
                                  'Color',MyBlue,...
                                  'LineWidth' , 3,...
                                  'MarkerSize', 25,...
                                  'MarkerFaceColor',MyRed,...
                                  'DisplayName', 'CE Optimum');


    set(gca,'FontName'   ,'Helvetica');
    set(gca,'FontSize'   ,18         );

    leg = [fig2(1); fig3; fig4; fig5];
    leg = legend(leg,'Location','SouthEast','FontSize',12);
    %leg = legend(leg,'Location','Best','FontSize',12);
    hold off;

    saveFileName = sprintf('ContourPlot_%s.eps', funcName);
    exportgraphics(gca, saveFileName, 'Resolution', 300);
end


% -----------------------------------------------------------------
% Interpretation of ExitFlag
% -----------------------------------------------------------------
function printExitFlagMeaning(ExitFlag)
    switch ExitFlag
        case 1
            disp(' Termination: Maximum number of iterations reached.');
        case 2
            disp(' Termination: Solution stalled.');
        case 3
            disp(' Termination: Maximum number of function evaluations reached.');
        case 4
            disp(' Termination: Function value variation is small .');
        case 5
            disp(' Termination: Standard deviation variation is small.');
        case 6
            disp(' Termination: Minimum function value criterion met.');
        otherwise
            disp(' Termination: Ended with an unknown exit flag.');
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% Benchmark functions
% -----------------------------------------------------------------

% Ackley function
function F = Ackley(x)
    n = size(x, 2);
    sum_sq = sum(x.^2, 2);
    sum_cos = sum(cos(2*pi*x), 2);
    F = -20 * exp(-0.2*sqrt(1/n * sum_sq)) - exp(1/n * sum_cos) + 20 + exp(1);
end

% Beale Function
function F = Beale(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = (1.5 - x1 + x1.*x2).^2 + (2.25 - x1 + x1.*x2.^2).^2 + (2.625 - x1 + x1.*x2.^3).^2;
end

% Booth function
function F = Booth(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = (x1 + 2*x2 - 7).^2 + (2*x1 + x2 - 5).^2;
end

% Bukin function N.6
function F = BukinN6(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = 100 * sqrt(abs(x2 - 0.01*x1.^2)) + 0.01*abs(x1 + 10);
end

% Cross-in-tray function
function F = CrossInTray(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = -0.0001 * (abs(sin(x1).*sin(x2).*exp(abs(100 - sqrt(x1.^2+x2.^2)/pi))) + 1).^0.1;
end

% Dixon-Price function
function F = DixonPrice(x)
    n = size(x, 2);
    term1 = (x(1) - 1)^2;
    sum_seq = (2:n)';
    terms = 2 * sum_seq .* ((2 * x(2:end).^2 - x(1:end-1)).^2);
    F = term1 + sum(terms, 2);
end

% Easom function
function F = Easom(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = -cos(x1).*cos(x2).*exp(-(x1-pi).^2-(x2-pi).^2);
end

% Eggholder function
function F = Eggholder(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = -(x2+47) .* sin(sqrt(abs(x2+x1/2+47))) - x1 .* sin(sqrt(abs(x1-(x2+47))));
end

% Goldstein-Price function
function F = GoldsteinPrice(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = (1 + (x1 + x2 + 1).^2 .* (19 - 14*x1 + 3*x1.^2 - 14*x2 + 6*x1.*x2 + 3*x2.^2)) .* ...
        (30 + (2*x1 - 3*x2).^2 .* (18 - 32*x1 + 12*x1.^2 + 48*x2 - 36*x1.*x2 + 27*x2.^2));
end

% Griewank function
function F = Griewank(x)
    [m, n] = size(x);
    i = repmat(1:n, m, 1);
    sum_sq = sum(x.^2 / 4000, 2);
    prod_cos = prod(cos(x ./ sqrt(i)), 2);
    F = sum_sq - prod_cos + 1;
end

% Himmelblau's function
function F = Himmelblau(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = (x1.^2 + x2 - 11).^2 + (x1 + x2.^2 - 7).^2;
end

% Holder Table Function
function F = HolderTable(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = -abs(sin(x1).*cos(x2).*exp(abs(1 - sqrt(x1.^2+x2.^2)/pi)));
end

% Levi function N.13
function F = LeviN13(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = sin(3*pi*x1).^2 + (x1-1).^2 .* (1+sin(3*pi*x2).^2) + (x2-1).^2 .* (1+sin(2*pi*x2).^2);
end

% Matyas function
function F = Matyas(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = 0.26*(x1.^2 + x2.^2) - 0.48*x1.*x2;
end

% McCormick function
function F = McCormick(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = sin(x1+x2) + (x1-x2).^2 - 1.5*x1 + 2.5*x2 + 1;
end

% Rastrigin function
function F = Rastrigin(x)
    n = size(x, 2);
    A = 10;
    F = A*n + sum(x.^2 - A * cos(2 * pi * x), 2);
end

% Rosenbrock function
function F = Rosenbrock(x)
    F = sum(100*(x(:,2:end) - x(:,1:end-1).^2).^2 + (1 - x(:,1:end-1)).^2, 2);
end

% Schaffer function N.2 
function F = SchafferN2(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = 0.5 + (sin(x1.^2 - x2.^2).^2 - 0.5) ./ (1 + 0.001*(x1.^2 + x2.^2)).^2;
end

% Schaffer function N.4
function F = SchafferN4(x)
    x1 = x(:,1);
    x2 = x(:,2);
    F = 0.5 + (cos(sin(abs(x1.^2 - x2.^2))).^2 - 0.5) ./ (1 + 0.001*(x1.^2 + x2.^2)).^2;
end

% Shekel function with random points
function F = Shekel(x)
    m = 10; % Number of local minima points
    A = [4 1 8 6 3 2 5 8 6 7; 4 1 8 6 3 2 5 8 6 7]';
    C = 1/m*[1 2 2 4 4 6 3 7 5 5]';
    F = zeros(size(x,1), 1);
    for i = 1:m
        F = F - 1 ./ (sum((x - A(i,:)).^2, 2) + C(i));
    end
end

% Sphere function
function F = Sphere(x)
    F = sum(x.^2, 2);
end

% Styblinski–Tang function
function F = StyblinskiTang(x)
    F = sum(x.^4 - 16*x.^2 + 5*x, 2) / 2;
end

% Three-hump Camel function
function F = ThreeHumpCamel(x)
    x1 = x(:,1);
    x2 = x(:,1);
    F = 2*x1.^2 - 1.05*x1.^4 + (x1.^6)/6 + x1.*x2 + x2.^2;
end

% Zakharov function
function F = Zakharov(x)
    n = size(x, 2);
    sum1 = sum(x.^2, 2);
    sum2 = sum(0.5 * (1:n) .* x, 2);
    F = sum1 + sum2.^2 + sum2.^4;
end
% -----------------------------------------------------------------