% -----------------------------------------------------------------
%  PlotTruss10.m
% -----------------------------------------------------------------
%  programmers: Marcos Vinicius Issa
%               Anderson Pereira
%               Americo Cunha Jr
%
%  Originally programmed in: Apr 04, 2024
%           Last updated in: Aug 29, 2024
% -----------------------------------------------------------------
%  Plot a 10 bars truss given the elements cross-section area.
% -----------------------------------------------------------------
function PlotTruss10(x,MyTruss,MyTitle)

    % truss structure parameters
    h     = MyTruss.h;
    NODES = MyTruss.NODES;
    ELEM  = MyTruss.ELEM;
    
    % custom color
    grayColor = [.7 .7 .7];
    
    figure('DefaultAxesFontSize',10)
    clf
    hold on
    
    % Support 1
    xl1 = [-1.0 -0.7]; 
    yl1 = [-1.0 1.0];
    [X1,Y1] = hatch_coordinates(xl1, yl1, 0.15) ;
    plot(X1,Y1,'k','linewidth',1.0);
    a11 = [0; -0.7; -0.7];
    a12 = [0; 1.0; -1.0];
    patch(a11,a12,grayColor,'EdgeColor','none');
    plot([a11(2) a11(3)],[a12(2) a12(3)],'k','linewidth',2.5);

    % Support 2
    xl2 = [-1.0 -0.7]; 
    yl2 = [h-1.0 h+1.0];
    [X2,Y2] = hatch_coordinates(xl2, yl2, 0.15) ;
    plot(X2,Y2,'k','linewidth',1.0);
    a21 = [0; -0.7; -0.7];
    a22 = [h; h+1; h-1];
    patch(a21,a22,grayColor,'EdgeColor','none');
    plot([a21(2) a21(3)],[a22(2) a22(3)],'k','linewidth',2.5);
   
    for i = 1:length(x)
        patch('Faces'          ,ELEM(i,:),...
              'Vertices'       ,NODES,...
              'EdgeColor'      ,grayColor,...
              'LineWidth'      ,x(i));

        patch('Faces'          ,ELEM(i,:),...
              'Vertices'       ,NODES,...
              'EdgeColor'      ,'none',...
              'FaceColor'      ,'none', ...
              'MarkerEdgeColor','blue',...
              'Marker'         ,'o',...
              'MarkerFaceColor','white',...
              'MarkerSize'     ,10,... 
              'LineWidth'      ,3);
    end
    
    for i = 1:4
        plot([NODES(i,1)],[NODES(i,2)]  , ...
            'Marker'                    , 'o', ...
            'MarkerEdgeColor'           ,[0.9290 0.6940 0.1250], ...
            'MarkerSize'                , 17, ... 
            'LineWidth'                 , 4.1);
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