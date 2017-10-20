function [H,S] = master_matrices_3D_Electronic(Mesh_info)

% Recover Mesh info
pattern_vec     = Mesh_info.pattern_vec;
Vs              = Mesh_info.Vs;
n               = Mesh_info.n;
elenodes        = Mesh_info.elenodes;
coordinates     = Mesh_info.coordinates;

n_eles = size(elenodes,1);


i_col = zeros((n^3)^2*n_eles,1);
i_row = zeros((n^3)^2*n_eles,1);
H_val = zeros((n^3)^2*n_eles,1);
S_val = zeros((n^3)^2*n_eles,1);

for i = 1:n_eles
    i
    dofs_in_ele = zeros(n^3,1);
    for j = 1:length(elenodes(i,:))
        dofs_in_ele(j) = elenodes(i,j);
    end
    
    % extract element potential
    V = Vs{pattern_vec(i)+1};
 
    [H_ele,S_ele] = ...
        element_matrices_3D_Electronic(n,coordinates(elenodes(i,:),:),V);
    
    % row and column indices for current element in master system
    i_row_ele = repmat(dofs_in_ele,[1,n^3]);
    i_col_ele = i_row_ele';
    
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    
    % master row and column index vectors
    i_col(((n^3)^2*(i-1)+1):((n^3)^2*i)) = i_col_ele;
    i_row(((n^3)^2*(i-1)+1):((n^3)^2*i)) = i_row_ele;
    
    % master Hamiltonian and Structure value vectors
    H_val(((n^3)^2*(i-1)+1):((n^3)^2*i)) = H_ele(:);
    S_val(((n^3)^2*(i-1)+1):((n^3)^2*i)) = S_ele(:);
    
end

H = sparse(i_row,i_col,H_val);
S = sparse(i_row,i_col,S_val);

H = (1/2)*(H+H');
S = (1/2)*(S+S');