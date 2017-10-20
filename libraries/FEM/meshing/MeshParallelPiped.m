function [coordinates,elenodes] = MeshParallelPiped(pattern,eleorder,R)


%% Create grid of points
n_ele_a1 = size(pattern,1);
n_ele_a2 = size(pattern,2);
n_ele_a3 = size(pattern,3);

n_node_a1 = n_ele_a1*eleorder + 1;
n_node_a2 = n_ele_a2*eleorder + 1;
n_node_a3 = n_ele_a3*eleorder + 1;

% number of nodes
n_nodes = n_node_a1*n_node_a2*n_node_a3;  

% number of dofs in periodic system
% n_dof_per(q1) = (n_node_a1-1)*(n_node_a2-1)*(n_node_a3-1);

% number of dofs before periodicity is enforced
% n_dof = n_nodes;

% number of elements
n_eles  = n_ele_a2*n_ele_a1*n_ele_a3;                    

% 3d grid in lattice vector coordinates
a1vec = linspace(0,1,n_node_a1);
a2vec = linspace(0,1,n_node_a2);
a3vec = linspace(0,1,n_node_a3);
[a1grid,a2grid,a3grid] = ndgrid(a1vec,a2vec,a3vec);

% vectorize lattice vector grid
a1locs = reshape(a1grid,n_nodes,1);
a2locs = reshape(a2grid,n_nodes,1);
a3locs = reshape(a3grid,n_nodes,1);

% % 3d grid in cartesian coordinates
% xgrid = (a1grid*a1(1) + a2grid*a2(1) + a3grid*a3(1));
% ygrid = (a1grid*a1(2) + a2grid*a2(2) + a3grid*a3(2));
% zgrid = (a1grid*a1(3) + a2grid*a2(3) + a3grid*a3(3));

% coordinates = [a1,a2,a3]*[a1locs',a2locs'
% R = [a1,a2,a3];
coordinates = [a1locs,a2locs,a3locs]*R';

% create node array
pattern_vec = reshape(pattern,n_eles,1);
node_vec = 1:n_nodes;
nodes = reshape(node_vec,n_node_a1,n_node_a2,n_node_a3);

% define element connectivity matrix
elenodes = zeros(n_eles,(eleorder+1)^3);
elenode_order = elenodes;

for i = 1:n_ele_a1
    for j = 1:n_ele_a2
        for k = 1:n_ele_a3
            x_bl = (i-1)*eleorder+1;
            y_bl = (j-1)*eleorder+1;
            z_bl = (k-1)*eleorder+1;

            elenode_cube = zeros((eleorder+1),(eleorder+1),(eleorder+1));
            for q = 1:(eleorder+1)
                elenode_cube(:,:,q) = nodes([x_bl:x_bl+eleorder],...
                    [y_bl:y_bl+eleorder],z_bl+q-1);
            end


            elenodes((i-1)*n_ele_a2*n_ele_a3+(j-1)*n_ele_a3+k,:) = ...
                elenode_cube(:);
        end
    end
end