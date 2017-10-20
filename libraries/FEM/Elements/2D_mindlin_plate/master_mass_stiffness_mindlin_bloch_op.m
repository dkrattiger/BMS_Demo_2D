function [K,M] = master_mass_stiffness_mindlin_bloch_op(Mesh_info)

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

field_names = {'s0','sx','sy','sxx','sxy','syy'};
for i = 1:length(field_names)
    K_val.(field_names{i}) = zeros((npn*n^2)^2*n_eles,1);
end

% K_valx = zeros((npn*n^2)^2*n_eles,1);
% K_valy = zeros((npn*n^2)^2*n_eles,1);
% K_valxy = zeros((npn*n^2)^2*n_eles,1);
% M_val = zeros((npn*n^2)^2*n_eles,1);

% K = zeros(n_dof);
% M = zeros(n_dof,n_dof);

for i = 1:n_eles
    i
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
        element_mass_stiffness_mindlin_bloch_op(Df,Dc,rho,h,n,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),elenode_order,eps);
    
%     K_ele.s0
%     K_ele.sx
%     K_ele.sxx
%     K_ele.sy
%     K_ele.syy
%     K_ele.sxy
%     pause

    format long g
%     A0 = K_ele.s0;
%     A1 = (K_ele.sx);
%     A2 = (K_ele.sxx);
%     L = sort(polyeig(A0,A1,A2))
%     figure(2014)
%     plot(sort(abs(L)))
%     pause
% A = ([A0,zeros(n_dof_per);zeros(n_dof_per),eye(n_dof_per)]);
% B = ([-A1,-A2;eye(n_dof_per),zeros(n_dof_per)]);
% 
% eigs(A,B,10,'sm')
%     eigs(K_ele.s0)
    
%     [K_ele_test,M_ele_test] = ...
%         element_mass_stiffness_mindlin(Df,Dc,rho,h,n,...
%         xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),elenode_order,eps);
%     eigs(K_ele_test)
%     pause
    % add element mass and stiffness into appropriate locations
%     M(dofs_in_ele,dofs_in_ele) = M(dofs_in_ele,dofs_in_ele) + M_ele;
%     K(dofs_in_ele,dofs_in_ele) = K(dofs_in_ele,dofs_in_ele) + K_ele;
    
%     col_indices = ones(n^2*npn,1)*dofs_in_ele;
%     row_indices = col_indices';
    
    row_indices = dofs_in_ele*ones(1,n^2*npn);
    col_indices = row_indices.';

    i_col(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = col_indices(:);
    i_row(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = row_indices(:);
%     K_val0(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = K_ele.s0(:);
%     K_valx(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = K_ele.sx(:);
%     K_valy(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = K_ele.sy(:);
%     K_valxy(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = K_ele.sxy(:);  
    for j = 1:length(field_names)
        K_val.(field_names{j})(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = ...
        K_ele.(field_names{j})(:);
    end
    M_val(((i-1)*(npn*n^2)^2+1):i*(npn*n^2)^2) = M_ele(:);
end

% form sparse matrices using row indices, column indices, and values
for j = 1:length(field_names)
    % form sparse matrix and symmetrize 
    K.(field_names{j}) = sparse(i_row,i_col,K_val.(field_names{j}));
%     K.(field_names{j}) = (1/2)*(K.(field_names{j})+K.(field_names{j})');
end
% eigs(K.s0,6,'sm')
% pause
% form sparse matrix then symmetrize
M = sparse(i_row,i_col,M_val);
M = (1/2)*(M+M');