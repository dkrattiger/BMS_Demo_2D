function K = StiffnessAssemble(X,ele_info,Dcell)

% ele_info should contain nodes in element, element material index, element
% type

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
zlocs           = Mesh_info.zlocs;

n_kap  = size(kappa,2);
n_eles = size(elenodes,1);

% K = zeros(n_dof,n_dof,n_kap);
% M = zeros(n_dof,n_dof);
% 
% K = spalloc(n_dof,n_dof,3*n^3*n_eles);
% M = spalloc(n_dof,n_dof,3*n^3*n_eles);

i_col = zeros(3*n^3*n_eles,1);
i_row = zeros(3*n^3*n_eles,1);
K_val = zeros(3*n^3*n_eles,1);
M_val = zeros(3*n^3*n_eles,1);

for i = 1:n_eles
    i
    dofs_in_ele = zeros(3*n^3,1);
    for j = 1:length(elenodes(i,:))
        en = elenodes(i,j);
        dofs_in_ele(3*(j-1)+1:3*j) = [3*(en-1)+1:3*en];
    end

    rho = rhos{pattern_vec(i)+1};
    D   = Ds{pattern_vec(i)+1};
%     pause
 
% elenodes(i,:)
% D,rho,n,kappa, xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),zlocs(elenodes(i,:)),...
%         elenode_order

    [K_ele,M_ele] = element_mass_stiffness_3DE(D,rho,n,kappa,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),zlocs(elenodes(i,:)),...
        elenode_order);
%     size(K_ele)
%     size(M_ele)
%     
%     figure(134)
%     subplot(1,2,1)
%     spy(K_ele)
%     subplot(1,2,2)
%     spy(M_ele)
    
%     lambdas = eigs(K_ele,M_ele,10,'sm');
%     if min(lambdas)<-1e-6*max(max(K_ele))
%         figure(4343);clf
%         plot(lambdas)
%         pause
%     end
    
    
    i_row_ele = repmat(dofs_in_ele,[1,3*n^3]);
    i_col_ele = i_row_ele';
    
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    i_col(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = i_col_ele;
    i_row(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = i_row_ele;
    
    K_val(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = K_ele(:);
    M_val(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = M_ele(:);
    
% %     dofs_in_ele = sort(dofs_in_ele);
%     % Master Mass
%     M(dofs_in_ele,dofs_in_ele) = M(dofs_in_ele,dofs_in_ele) + M_ele;
%     
%     % Master Stiffness
% %     for j = 1:n_kap
%     K(dofs_in_ele,dofs_in_ele,:) = K(dofs_in_ele,dofs_in_ele,:)+...
%         K_ele;
% %     end
%     
end

K = sparse(i_row,i_col,K_val);
M = sparse(i_row,i_col,M_val);