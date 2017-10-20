function [K,M] = master_mass_stiffness_planeStrainQuad_BlochOperator(kappa,Mesh_info)

% Recover Mesh info
Ds              = Mesh_info.Ds;
rhos            = Mesh_info.rhos;
pattern_vec     = Mesh_info.pattern_vec;
n               = Mesh_info.n;
n_dof           = Mesh_info.n_dof;
elenodes        = Mesh_info.elenodes;
elenode_order   = Mesh_info.elenode_order;
xlocs           = Mesh_info.xlocs;
ylocs           = Mesh_info.ylocs;

n_kap  = size(kappa,2);
n_eles = size(elenodes,1);

% K = zeros(n_dof);
% M = zeros(n_dof,n_dof);
% i_col_vec = zeros(n_eles*
i_col = zeros(2*n^2*n_eles,1);
i_row = zeros(2*n^2*n_eles,1);
K_val = zeros(2*n^2*n_eles,1);
M_val = zeros(2*n^2*n_eles,1);
if n_kap>1
    K_valx = zeros(2*n^2*n_eles,1);
    K_valy = zeros(2*n^2*n_eles,1);
    K_valxx = zeros(2*n^2*n_eles,1);
    K_valyy = zeros(2*n^2*n_eles,1);
    K_valxy = zeros(2*n^2*n_eles,1);
end

for i = 1:n_eles
    i
    dofs_in_ele = zeros(2*n^2,1);
    for j = 1:length(elenodes(i,:))
        en = elenodes(i,j);
        dofs_in_ele(2*j-1:2*j) = [2*en-1,2*en];
    end

    rho = rhos{pattern_vec(i)+1};
    D   = Ds{pattern_vec(i)+1};
    
    [K_ele,Kx_ele,Ky_ele,Kxx_ele,Kyy_ele,Kxy_ele,M_ele] = ...
        element_mass_stiffness_planeStrainQuad_BlochOperator(D,rho,n,kappa,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),elenode_order);
    
    i_row_ele = repmat(dofs_in_ele,[1,2*n^2]);
    i_col_ele = i_row_ele';
    
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    i_col(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = i_col_ele;
    i_row(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = i_row_ele;
    
    K_val(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = K_ele(:);
    M_val(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = M_ele(:);
    
    % if Bloch operator
    if n_kap>1
        K_valx(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = Kx_ele(:);
        K_valy(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = Ky_ele(:);
        K_valxx(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = Kxx_ele(:);
        K_valxy(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = Kxy_ele(:);
        K_valyy(((2*n^2)^2*(i-1)+1):((2*n^2)^2*i)) = Kyy_ele(:);
    end
end


M = sparse(i_row,i_col,M_val);
M = (1/2)*(M+M');

if n_kap>1
    K.s0 = sparse(i_row,i_col,K_val);
    K.sx = sparse(i_row,i_col,K_valx);
    K.sy = sparse(i_row,i_col,K_valy);
    K.sxx = sparse(i_row,i_col,K_valxx);
    K.sxy = sparse(i_row,i_col,K_valxy);
    K.syy = sparse(i_row,i_col,K_valyy);
    
    K.s0 = (1/2)*(K.s0+K.s0');
    K.sx = (1/2)*(K.sx+K.sx');
    K.sy = (1/2)*(K.sy+K.sy');
    K.sxx = (1/2)*(K.sxx+K.sxx');
    K.sxy = (1/2)*(K.sxy+K.sxy');
    K.syy = (1/2)*(K.syy+K.syy');
else
    K = sparse(i_row,i_col,K_val);
    K = (1/2)*(K+K');
end

