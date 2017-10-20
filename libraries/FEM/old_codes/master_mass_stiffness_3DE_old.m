function [K,M] = master_mass_stiffness_3DE(kappa,Mesh_info)

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

K = zeros(n_dof,n_dof,n_kap);
M = zeros(n_dof,n_dof);

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
 
    [K_ele,M_ele] = element_mass_stiffness_3DE(D,rho,n,kappa,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),zlocs(elenodes(i,:)),...
        elenode_order);
    
%     dofs_in_ele = sort(dofs_in_ele);
    % Master Mass
    M(dofs_in_ele,dofs_in_ele) = M(dofs_in_ele,dofs_in_ele) + M_ele;
    
    % Master Stiffness
%     for j = 1:n_kap
    K(dofs_in_ele,dofs_in_ele,:) = K(dofs_in_ele,dofs_in_ele,:)+...
        K_ele;
%     end
    
end