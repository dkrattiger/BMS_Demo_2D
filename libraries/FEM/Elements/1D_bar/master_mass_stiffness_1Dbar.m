function [K,M] = master_mass_stiffness(pattern,E12,rho12,L,A,kappa,n)
    rhos   = rho12(pattern+1);
    Es = E12(pattern+1);

    % unit cell model properties
    n_eles  = length(pattern);
    n_nodes = (n-1)*n_eles+1;
    n_dof = n_nodes;
    n_kap   = length(kappa);

    % node locations (x coordinate only)
    node_locs = linspace(0,L,n_nodes);

    % create map that links the 1st and last node, so that the system 
    % will behave periodically
%     node_map = [1:(n_nodes-1),1];

    % nodes in each element
    ele_nodes = zeros(n_eles,n);
    for i = 1:n
        ele_nodes(:,i) = linspace(i,i+(n-1)*(n_eles-1),n_eles)';
    end

    % preallocate mass and stiffness arrays
    M = zeros(n_dof,n_dof);
    K = zeros(n_dof,n_dof,n_kap);
    M_ele = zeros(n,n,n_eles);
    K_ele = zeros(n,n,n_kap,n_eles);

    % loop through elements and integrate shape functions, then create 
    % element and master stiffnesses and masses

    for i = 1:n_eles

        % Find Element Stiffness
        [K_ele(:,:,:,i),M_ele(:,:,i)] = element_mass_stiffness(...
            Es(i),A,rhos(i),n,kappa,node_locs(ele_nodes(i,:)));

        % Form master mass
%         index = node_map(ele_nodes(i,:));
        index = ele_nodes(i,:);
        M(index,index) = M(index,index) + M_ele(:,:,i);


        % Form master stiffness
        for j = 1:n_kap
            % form master stiffness
            K(index,index,j) = K(index,index,j) + K_ele(:,:,j,i);
        end

    end