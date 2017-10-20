function [K_ele,M_ele] = element_mass_stiffness_3Dto1D_SAFE(D,rho,n,X,elenode_order)

    n_nodes = n;
    ndof = 3*n_nodes;
    
    % Function lgwt generates gauss legendre points and weights
    [gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);

    
    % Preallocate arrays
    M_ele = zeros(ndof,ndof);
    K_ele.s0 = zeros(ndof,ndof);
    K_ele.sy = zeros(ndof,ndof);
    K_ele.sz = zeros(ndof,ndof);
    K_ele.syy = zeros(ndof,ndof);
    K_ele.szz = zeros(ndof,ndof);
    K_ele.syz = zeros(ndof,ndof);
    
    % node locations in isoparametric coordinates
    zetas = linspace(-1,1,n);
    
    % loop over gauss quadriture points to perform integration of shape 
    % functions, then create element stiffness and mass
    jacobsum = 0;
    for i = 1:length(gq_pts) 

        zeta = gq_pts(i);
        weightz = gq_wts(i);  

        N = zeros(1,n_nodes);
        dNdz = zeros(1,n_nodes);
        
        iis = [1:n];
        for k = 1:n_nodes
            q = elenode_order(k);


            % index for mesh
            ii = iis(q);

            % zeta polynomial
            rz = zetas; rz(ii) = [];
            pz = poly(rz);

            % calculate normalization coefficient
            pnorm = polyval(pz,zetas(ii));

            % shape function
            N(k) = polyval(pz,zeta)/pnorm;

            % shape function derivatives
            dpdz = pz(1:end-1).*((length(pz)-1):-1:1);
            dNdz(k) = polyval(dpdz,zeta)/pnorm;
        end

        % collect shape function derivatives
        dNdZ = dNdz;

        J = [dNdz*X(:,1)];

        % compute real space derivatives
        dNdX = J\dNdZ;

%            
        jacob = det(J);
        jacobsum = jacobsum+jacob;

        % Shape Function Matrix
        Nmat = zeros(3,ndof);
        for k = 1:n_nodes
            Nmat(:,(3*(k-1)+1:3*k)) = N(k)*eye(3);
        end

        % Element Mass Matrix
        m_add = weightz*jacob*rho*(Nmat.'*Nmat);
        M_ele = M_ele+m_add;      

        % strain-displacement derivative locations
        Lx = [1,    0,  0;...
              0,    0,  0;...
              0,    0,  0;...
              0,    0,  0;...
              0,    0,  1;...
              0,    1,  0];

        Ly = [0,    0,  0;...
              0,    1,  0;...
              0,    0,  0;...
              0,    0,  1;...
              0,    0,  0;...
              1,    0,  0];

        Lz = [0,    0,  0;...
              0,    0,  0;...
              0,    0,  1;...
              0,    1,  0;...
              1,    0,  0;...
              0,    0,  0];


        % Element Stiffness Matrix
        Bx = zeros(6,ndof);
        By = zeros(6,ndof);
        Bz = zeros(6,ndof);

        for k1 = 1:n_nodes
            Bx(:,(3*k1-2:3*k1)) = Lx*dNdX(1,k1);
            By(:,(3*k1-2:3*k1)) = Ly*N(1,k1);
            Bz(:,(3*k1-2:3*k1)) = Lz*N(1,k1);
        end
        
        % Contribution to stiffness terms from current GQ point
        k_add0 = weightz*jacob*(Bx'*D*Bx);
        k_addy = 1i*weightz*jacob*(Bx'*D*By-By'*D*Bx);
        k_addz = 1i*weightz*jacob*(Bx'*D*Bz-Bz'*D*Bx);
        k_addyy = weightz*jacob*(By'*D*By);
        k_addzz = weightz*jacob*(Bz'*D*Bz);
        k_addyz = weightz*jacob*(Bz'*D*By+By'*D*Bz);

        % sum (integrate) contribution from current GQ point
        K_ele.s0 = K_ele.s0+k_add0;
        K_ele.sy = K_ele.sy+k_addy;
        K_ele.sz = K_ele.sz+k_addz;
        K_ele.syy = K_ele.syy+k_addyy;
        K_ele.szz = K_ele.szz+k_addzz;
        K_ele.syz = K_ele.syz+k_addyz;
    end
    
    % check jacobian
    if jacobsum<0
        disp('negative jacobian')
        K_ele.s0 = -K_ele.s0;
        K_ele.sy = -K_ele.sy;
        K_ele.sz = -K_ele.sz;
        K_ele.syy = -K_ele.syy;
        K_ele.szz = -K_ele.szz;
        K_ele.syz = -K_ele.syz;
        M_ele = -M_ele;
    end