function [h_patch,h_line] = plot_FEM_model(coordinates,patchfaces,C,edges)

%
% Dimitri Krattiger
%
% Description
% ===========
% This function plots a finite element model whose elements are specified
% with "patchfaces" and whose nodes are specified in "coordinates"
%
% inputs
% ======
% coordinates   = node coordinates
%
% patchfaces    = array of node indices that define element patches
%
% C             = array of colors specifying patch colors
%
% edges         = array specifying lines that mark model feature edges
%
% outputs
% =======
% h_patch       = handle to patches
% 
% h_line        = handle to feature edge lines


% number of dimensions
[n_nodes,n_dim] = size(coordinates);

% plot dots for nodes?
plot_nodes = false;
plot_corner_nodes = false;
n = size(patchfaces,2)/4;
if plot_corner_nodes
    i_node_plot = patchfaces(:,1:n:end);
else
    i_node_plot = patchfaces(:);
end
i_node_plot = unique(i_node_plot(:));
i_node_plot(isnan(i_node_plot)) = [];
   

% plot patches for each finite element in unit cell
h_patch = patch('faces',patchfaces,...
    'vertices',coordinates,'facecolor','w');
set(h_patch,'edgecolor',0.85*[1,1,1])

% define patch coloring
if ~isempty(C)
    %C = C/max(abs(C(:)));
    set(h_patch,'FaceColor','flat')
    set(h_patch,'FaceVertexCData',C);
end

hold on

if n_dim == 2
    
    x = coordinates(:,1);
    y = coordinates(:,2); 
    
    % plot feature edges
    h_line = [];    
    for i = 1:length(edges)
        h_line(i) = plot(x(edges{i})',y(edges{i})');
    end
    
    % define axis labels and scaling
    xlabel('x (m)');ylabel('y (m)');
    
    if plot_nodes
        plot(x(i_node_plot)',y(i_node_plot)',...
            'k.','markersize',10)
    end
    
elseif n_dim == 3
    
    
    x = coordinates(:,1);
    y = coordinates(:,2);
    z = coordinates(:,3);
    
    % plot feature edges
    h_line = [];
    if iscell(edges)
        for i = 1:length(edges)
            h_line(i) = plot3(x(edges{i})',y(edges{i})',z(edges{i})');
        end
    else
        h_line = plot3(x(edges)',y(edges)',z(edges)');
    end
    
    % define axis labels and scaling
    xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');
    
    if plot_nodes
        plot3(x(i_node_plot)',y(i_node_plot)',z(i_node_plot)',...
            'k.','markersize',10)
    end
end
    
% format model edges
for i = 1:length(edges)
    set(h_line(i),'color','k')
    set(h_line(i),'linestyle','-')
    set(h_line(i),'linewidth',2)
end

% specify renderer
set(gcf,'renderer','opengl')
hold off

% % use persepective visualization
% camproj('perspective')
% campos([-Lx,-Ly,Lz*1.2]*0.7);