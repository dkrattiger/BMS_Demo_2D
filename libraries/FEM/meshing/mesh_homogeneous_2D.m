function [xlocs,ylocs,elenodes,pattern,patchfaces,C] = ...
     mesh_homogeneous_2D(n,n_ele_side,L)
         

n_nodes_side = n_ele_side*(n-1)+1;
% n_nodes_incl = n_ele_incl*(n-1)+1;
n_ele = n_ele_side^2;
n_nodes = n_nodes_side^2;

%% cube nodes;
% ======================================================================= %
vec = linspace(0,L,n_nodes_side);
[xgrid,ygrid] = meshgrid(vec,vec);

xlocs = xgrid(:);
ylocs = ygrid(:);


%% define element connectivity matrix
% ======================================================================= %
elenodes = zeros(n_ele,n^2);
nodes_square = reshape(1:n_nodes,n_nodes_side,n_nodes_side);

corners = [1,n,n^2,n^2-n+1];
sides = [2:n-1,2*n:n:n^2-n,(n^2-1:-1:n^2-n+2),(n^2-2*n+1:-n:n+1)];
interior = 1:n^2; interior([corners,sides]) = [];
elenode_order = [corners,sides,interior];

for i = 1:n_ele_side
    for j = 1:n_ele_side
        
        y_bl = (i-1)*(n-1)+1;
        x_bl = (j-1)*(n-1)+1;

        elenodes_temp = nodes_square([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)])';
        elenodes((i-1)*n_ele_side+j,:) = ...
            elenodes_temp(elenode_order);
    end
end

%% Define pattern (identifying element materials)
% ======================================================================= %
pattern = ones(n_ele_side,n_ele_side);
pattern = pattern(:);

%% Create patches for element faces and plot
% ======================================================================= %
ele_patch_edge = [1,5+0*(n-2):4+1*(n-2),...
                  2,5+1*(n-2):4+2*(n-2),...
                  3,5+2*(n-2):4+3*(n-2),...
                  4,5+3*(n-2):4+4*(n-2)];
ele_patch_full = [ele_patch_edge,4*(n-1)+1:n^2];
      
% six_patches = [bottom;top;side12;side23;side34;side41];
patchfaces = zeros(n_ele,4*(n-1));
patchfaces2 = zeros(n_ele,n^2);
[~,elenode_order_reverse] = sort(elenode_order);
for i = 1:n_ele
    patchfaces(i,:) = elenodes(i,ele_patch_edge);
    patchfaces2(i,:) = elenodes(i,elenode_order_reverse);
end

% coloring vector
C = ([pattern,pattern,pattern]+1)/3;

% node coordinate array
xylocs = [xlocs,ylocs];
