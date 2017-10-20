function [K,M] = master_mass_stiffness_Tet10(Mesh_info)

% Recover Mesh info
Ds              = Mesh_info.Ds;
rhos            = Mesh_info.rhos;
pattern_vec     = Mesh_info.pattern_vec;
% n               = Mesh_info.n;
% n_dof           = Mesh_info.n_dof;
elenodes        = Mesh_info.elenodes;
% elenode_order   = Mesh_info.elenode_order;
X               = Mesh_info.X;
% xlocs           = Mesh_info.xlocs;
% ylocs           = Mesh_info.ylocs;
% zlocs           = Mesh_info.zlocs;

% X = [xlocs,ylocs,zlocs];
% n_kap  = size(kappa,2);
n_eles = size(elenodes,1);

% K = zeros(n_dof,n_dof,n_kap);
% M = zeros(n_dof,n_dof);
% 
% K = spalloc(n_dof,n_dof,3*n^3*n_eles);
% M = spalloc(n_dof,n_dof,3*n^3*n_eles);

i_col = zeros(30^2*n_eles,1);
i_row = zeros(30^2*n_eles,1);
K_val = zeros(30^2*n_eles,1);
M_val = zeros(30^2*n_eles,1);

for i = 1:n_eles
    i
    dofs_in_ele = zeros(30,1);
    for j = 1:10
        en = elenodes(i,j);
        dofs_in_ele(3*(j-1)+1:3*j) = [3*(en-1)+1:3*en];
    end

    rho = rhos{pattern_vec(i)};
    D   = Ds{pattern_vec(i)};
    K_ele = Tet10Stiffness(X(elenodes(i,:),:),D);
    M_ele = Tet10Mass(X(elenodes(i,:),:),rho);    
    
    i_row_ele = repmat(dofs_in_ele,[1,30]);
    i_col_ele = i_row_ele';
    
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    i_col(((30)^2*(i-1)+1):((30)^2*i)) = i_col_ele;
    i_row(((30)^2*(i-1)+1):((30)^2*i)) = i_row_ele;
    
    K_val(((30)^2*(i-1)+1):((30)^2*i)) = K_ele(:);
    M_val(((30)^2*(i-1)+1):((30)^2*i)) = M_ele(:);

end

K = sparse(i_row,i_col,K_val);
M = sparse(i_row,i_col,M_val);