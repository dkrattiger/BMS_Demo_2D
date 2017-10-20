function [H0,Hx,Hy,Hz,S] = master_matrices_3D_Electronic_Bloch_Operator(kappa,Mesh_info)

% Recover Mesh info
pattern_vec     = Mesh_info.pattern_vec;
Vs              = Mesh_info.Vs;
n               = Mesh_info.n;
n_dof           = Mesh_info.n_dof;
elenodes        = Mesh_info.elenodes;
elenode_order   = Mesh_info.elenode_order;
xlocs           = Mesh_info.xlocs;
ylocs           = Mesh_info.ylocs;
zlocs           = Mesh_info.zlocs;
a0              = Mesh_info.a0;

n_kap  = size(kappa,2);
n_eles = size(elenodes,1);

H0 = zeros(n_dof,n_dof);
Hx = zeros(n_dof,n_dof,n_kap);
Hy = zeros(n_dof,n_dof,n_kap);
Hz = zeros(n_dof,n_dof,n_kap);
S = zeros(n_dof,n_dof);

for i = 1:n_eles
    i
    dofs_in_ele = zeros(n^3,1);
    for j = 1:length(elenodes(i,:))
        dofs_in_ele(j) = elenodes(i,j);
    end
    
    % extract element potential
    V = Vs{pattern_vec(i)+1};
 
    [H_ele_0,H_ele_kx,H_ele_ky,H_ele_kz,S_ele] = ...
        element_matrices_3D_Electronic_Bloch_Operator(n,kappa,...
        xlocs(elenodes(i,:)),ylocs(elenodes(i,:)),zlocs(elenodes(i,:)),...
        V,elenode_order);
    
%     dofs_in_ele = sort(dofs_in_ele);
    % Master Mass
    S(dofs_in_ele,dofs_in_ele) = S(dofs_in_ele,dofs_in_ele) + S_ele;
    
    % Master Stiffness
    H0(dofs_in_ele,dofs_in_ele,:) = H0(dofs_in_ele,dofs_in_ele,:)+...
        H_ele_0;
    Hx(dofs_in_ele,dofs_in_ele,:) = Hx(dofs_in_ele,dofs_in_ele,:)+...
        H_ele_kx;
    Hy(dofs_in_ele,dofs_in_ele,:) = Hy(dofs_in_ele,dofs_in_ele,:)+...
        H_ele_ky;
    Hz(dofs_in_ele,dofs_in_ele,:) = Hz(dofs_in_ele,dofs_in_ele,:)+...
        H_ele_kz;
%     Hv(dofs_in_ele,dofs_in_ele,:) = Hv(dofs_in_ele,dofs_in_ele,:)+...
%         H_ele_V;
end