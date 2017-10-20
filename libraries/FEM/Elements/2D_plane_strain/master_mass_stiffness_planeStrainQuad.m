function [K,M] = master_mass_stiffness_planeStrainQuad(Mesh_info)

% Recover Mesh info from structure
Ds              = Mesh_info.Ds;
rhos            = Mesh_info.rhos;
pattern_vec     = Mesh_info.pattern_vec;
n               = Mesh_info.n;
elenodes        = Mesh_info.elenodes;
xlocs           = Mesh_info.xlocs;
ylocs           = Mesh_info.ylocs;

% number of elements
n_eles = size(elenodes,1);

% preallocate vectors for sparse mass and stiffness matrix assembly
i_col = zeros(2*n^2*n_eles,1);
i_row = zeros(2*n^2*n_eles,1);
K_val = zeros(2*n^2*n_eles,1);
M_val = zeros(2*n^2*n_eles,1);

% loop through elements
for i = 1:n_eles
    
    % get element DOF indices
    dofs_in_ele = zeros(2*n^2,1);
    for j = 1:length(elenodes(i,:))
        en = elenodes(i,j);
        dofs_in_ele(2*j-1:2*j) = [2*en-1,2*en];
    end

    % element density
    rho = rhos{pattern_vec(i)+1};
    
    % element elasticity matrix (3x3)
    D   = Ds{pattern_vec(i)+1};
    
    % integrate element mass and stiffness matrices
    [K_ele,M_ele] = ...
        element_mass_stiffness_planeStrainQuad(D,rho,n,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)));
    
    % determine global indices of element DOFs
    i_row_ele = repmat(dofs_in_ele,[1,2*n^2]);
    i_col_ele = i_row_ele';
    
    % vectorize
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    % place element global indices into global index list
    i_col(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = i_col_ele;
    i_row(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = i_row_ele;
    
    % place element stiffness and mass matrix values into global mass and
    % stiffness matrix lists
    K_val(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = K_ele(:);
    M_val(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = M_ele(:);
    
end

% form sparse mass and stiffness matrices and ensure symmetry
M = sparse(i_row,i_col,M_val);
M = (1/2)*(M+M');
K = sparse(i_row,i_col,K_val);
K = (1/2)*(K+K');