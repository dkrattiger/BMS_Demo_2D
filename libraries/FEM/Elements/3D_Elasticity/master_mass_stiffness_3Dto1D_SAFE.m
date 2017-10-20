function [K,M] = master_mass_stiffness_3Dto1D_SAFE(Mesh_info)

% Recover Mesh info
Ds              = Mesh_info.Ds;
rhos            = Mesh_info.rhos;
pattern_vec     = Mesh_info.pattern_vec;
n               = Mesh_info.n;
elenodes        = Mesh_info.elenodes;
elenode_order   = Mesh_info.elenode_order;
X               = Mesh_info.X;

% number of elements
n_eles = size(elenodes,1);

% stiffness matrix fieldnames
fieldnames_K = {'s0','sy','sz','syy','szz','syz'};

i_col = zeros((3*n)^2*n_eles,1);
i_row = zeros((3*n)^2*n_eles,1);
M_val = zeros((3*n)^2*n_eles,1);
for j = 1:length(fieldnames_K)
    K_val.(fieldnames_K{j}) = zeros((3*n)^2*n_eles,1);
end
for i = 1:n_eles
    i
    dofs_in_ele = zeros(3*n,1);
    for j = 1:length(elenodes(i,:))
        en = elenodes(i,j);
        dofs_in_ele(3*(j-1)+1:3*j) = [3*(en-1)+1:3*en];
    end

    rho = rhos{pattern_vec(i)};
    D   = Ds{pattern_vec(i)};

    [K_ele,M_ele] = element_mass_stiffness_3Dto1D_SAFE(D,rho,n,...
                    X(elenodes(i,:),:),elenode_order);
    
    
    i_row_ele = repmat(dofs_in_ele,[1,3*n]);
    i_col_ele = i_row_ele';
    
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    i_col(((3*n)^2*(i-1)+1):((3*n)^2*i)) = i_col_ele;
    i_row(((3*n)^2*(i-1)+1):((3*n)^2*i)) = i_row_ele;
    
    
    M_val(((3*n)^2*(i-1)+1):((3*n)^2*i)) = M_ele(:);
    
    for j = 1:length(fieldnames_K)
        K_val.(fieldnames_K{j})(((3*n)^2*(i-1)+1):((3*n)^2*i)) = K_ele.(fieldnames_K{j})(:);
    end
end

% form sparse matrices and ensure
for j = 1:length(fieldnames_K)
    K.(fieldnames_K{j}) = sparse(i_row,i_col,K_val.(fieldnames_K{j}));
end
M = sparse(i_row,i_col,M_val);

% (NOT SURE IF THESE SHOULD BE SYMMETRIC SO DONT SYMMETRIZE)
% K = (1/2)*(K+K');
% M = (1/2)*(M+M');