function [K,M] = master_mass_stiffness_mindlin(Mesh_info)

verbose = false;

% Recover Mesh info
Dfs             = Mesh_info.Dfs;
Dcs             = Mesh_info.Dcs;
rhos            = Mesh_info.rhos;
hs              = Mesh_info.hs;
pattern_vec     = Mesh_info.pattern_vec;
n               = Mesh_info.n;
n_dof           = Mesh_info.n_dof;
elenodes        = Mesh_info.elenodes;
elenode_order   = Mesh_info.elenode_order;
xlocs           = Mesh_info.xlocs;
ylocs           = Mesh_info.ylocs;
epss            = Mesh_info.epss;

n_eles = size(elenodes,1); 
npn = 3; %number of dofs per node

i_col = zeros((npn*n^2)^2*n_eles,1);
i_row = zeros((npn*n^2)^2*n_eles,1);
K_val = zeros((npn*n^2)^2*n_eles,1);
M_val = zeros((npn*n^2)^2*n_eles,1);

% K = zeros(n_dof);
% M = zeros(n_dof,n_dof);

for i = 1:n_eles
    if verbose
        i
    end
    dofs_in_ele = zeros(npn*n^2,1);
    for j = 1:length(elenodes(i,:))
        en = elenodes(i,j);
        dofs_in_ele((npn*(j-1)+1):npn*j) = (npn*(en-1)+1):npn*en;
    end

    rho = rhos{pattern_vec(i)+1};
    Df   = Dfs{pattern_vec(i)+1};
    Dc   = Dcs{pattern_vec(i)+1};
    h    = hs{pattern_vec(i)+1};
    eps  = epss{pattern_vec(i)+1};
    
    % form element mass and stiffness
    [K_ele,M_ele] = ...
        element_mass_stiffness_mindlin(Df,Dc,rho,h,n,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),elenode_order,eps);

    % add element mass and stiffness into appropriate locations
%     M(dofs_in_ele,dofs_in_ele) = M(dofs_in_ele,dofs_in_ele) + M_ele;
%     K(dofs_in_ele,dofs_in_ele) = K(dofs_in_ele,dofs_in_ele) + K_ele;
    
%     col_indices = ones(n^2*npn,1)*dofs_in_ele;
%     row_indices = col_indices';
    
    row_indices = dofs_in_ele*ones(1,n^2*npn);
    col_indices = row_indices.';

    i_col(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = col_indices(:);
    i_row(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = row_indices(:);
    K_val(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = K_ele(:);  
    M_val(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = M_ele(:);
end

% form sparse matrices using row indices, column indices, and values
K = sparse(i_row,i_col,K_val);
M = sparse(i_row,i_col,M_val);

% symmetrize matrices
K = (1/2)*(K+K');
M = (1/2)*(M+M');