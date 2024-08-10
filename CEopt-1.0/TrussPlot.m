% -----------------------------------------------------------------
%  TrussPlot.m
% -----------------------------------------------------------------
%  programmer: Marcos Vinicius Issa
%              marcos.issa@uerj.br
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Jul 31, 2024
% -----------------------------------------------------------------
%  ﻿Plot a truss structure given the elements cross-section area.
% -----------------------------------------------------------------
function TrussPlot(x,MyTitle)

    % inch to meter convertion factor
    Inch2Meter = 0.0254;
    
    % truss geometric dimensions (m)
    l1 = 360*Inch2Meter; % 1st length
    l2 = 360*Inch2Meter; % 2nd length
    h  = 360*Inch2Meter; % heigth
    
    % mesh points in global coordinates (Nnodes x 2)
    XY_Mesh = [ (l1+l2)  h;
                 (l1+l2) 0;
                     l1  h;
                      l1 0;
                      0 l1;
                      0 0];
    
    % global coordinates of element nodes (Nelem x Nlnodes)
    GlobalNodeID = [3 5;
                    1 3;
                    4 6;
                    2 4;
                    3 4;
                    1 2;
                    4 5;
                    3 6;
                    2 3;
                    1 4];
    
    grayColor = [.7 .7 .7];
    
    figure('DefaultAxesFontSize',10)
    clf
    hold on
    
    % Support 1
    xl1 = [-1.0 -0.7] ; yl1 = [-1.0 1.0];
    [X1,Y1] = hatch_coordinates(xl1, yl1, 0.15) ;
    plot(X1,Y1,'red','linewidth',1.0);
    a11 = [0; -0.7; -0.7];
    a12 = [0; 1.0; -1.0];
    patch(a11,a12,grayColor,'EdgeColor','none');
    plot([a11(2) a11(3)],[a12(2) a12(3)],'red','linewidth',2.5);

    % Support 2
    xl2 = [-1.0 -0.7] ; yl2 = [h-1.0 h+1.0];
    [X2,Y2] = hatch_coordinates(xl2, yl2, 0.15) ;
    plot(X2,Y2,'red','linewidth',1.0);
    a21 = [0; -0.7; -0.7];
    a22 = [h; h+1; h-1];
    patch(a21,a22,grayColor,'EdgeColor','none');
    plot([a21(2) a21(3)],[a22(2) a22(3)],'red','linewidth',2.5);

    for i = 1:length(x)
        patch('Faces'          ,GlobalNodeID(i,:),...
              'Vertices'       ,XY_Mesh,...
              'EdgeColor'      ,grayColor,...
              'FaceColor'      ,'blue', ...
              'LineWidth'      ,x(i));

        patch('Faces'          ,GlobalNodeID(i,:),...
              'Vertices'       ,XY_Mesh,...
              'EdgeColor'      ,'none',...
              'FaceColor'      ,'none', ...
              'MarkerEdgeColor','blue',...
              'Marker'         ,'o',...
              'MarkerFaceColor','white',...
              'MarkerSize'     ,10,... 
              'LineWidth'      ,3);
    end


    set(gca,'xtick',[],'ytick',[])
    set(gca,'XColor', 'none','YColor','none')
    title(MyTitle,'FontSize',18)
    pause(1)
end
%------------------------------------------------------------

%------------------------------------------------------------
function [X,Y] = hatch_coordinates(xlim, ylim, xstep, ystep, merge)
    %// function [X,Y] = hatch_coordinates(xlim,ylim,xstep,ystep,merge)
    %//
    %// Return coordinates for plotting a hatch pattern
    %// The angle of the lines can be adjusted by varying the 
    % ratio xstep/ystep
    %// xlim and ylim are vectors with two elements, 
    % where the first element needs to be smaller than the second.

    % // set default options
    if nargin < 3 ; xstep = 1     ; end
    if nargin < 4 ; ystep = xstep ; end
    if nargin < 5 ; merge = true  ; end

    % // define base grid
    xpos = xlim(1):xstep:xlim(2) ; nx = numel(xpos) ;
    ypos = ylim(1):ystep:ylim(2) ; ny = numel(ypos) ;

    % // Create the coordinates
    nanline = NaN*ones(1,nx+ny-3) ;
    X = [ [ xpos(1)*ones(1,ny-2) xpos(1:end-1) ] ; ...
      [ xpos(2:end) xpos(end)*ones(1,ny-2) ] ; ...
      nanline ] ;
    Y = [ [ypos(end-1:-1:1) ylim(1)*ones(1,nx-2) ]  ; ...
      [ypos(end)*ones(1,nx-1) ypos(end-1:-1:2)] ; ...
      nanline ] ;

    % // merge if asked too
    if merge
        X = X(:) ;
        Y = Y(:) ;
    end
end
%------------------------------------------------------------