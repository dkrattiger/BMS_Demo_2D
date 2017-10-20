function [K_ele,M_ele] ...
    = element_mass_stiffness_mindlin(Df,Dc,rho,h,n,...
    xlocs,ylocs,elenode_order,eps)

    alpha = 5/6; % shear parameter
    
    % mass weighting matrix
    rho_mat = rho*diag([h,h^3/12,h^3/12]);
    %     rho_mat = rho*diag([h^3/12,h^3/12,h]);

    
    % element type
    npn = 3; % number of degrees of freedom per node
    n_nodes = n^2;  % number of nodes per element
    ndof = npn*n_nodes;   % number of degrees of freedom per element
    
    % Function lgwt generates gauss legendre points and weights
    [gq_pts,gq_wts] = lgwt(n,-1,1);
    
    
    % Preallocate arrays
%     K_ele = zeros(ndof);
    KB_ele = zeros(ndof);
    KS_ele_full = zeros(ndof);
    KS_ele_red = zeros(ndof);
    M_ele = zeros(ndof);

    % loop over gauss quadriture points to perform integration of shape 
    % functions, then create element stiffness and mass
    
    zetas = linspace(-1,1,n);
    etas = linspace(-1,1,n);
    
    %% Perform Full Integration of Shear Stiffness, Bending Stiffness, and Mass
    jacobsum = 0;
    for j = 1:length(gq_pts)
        for i = 1:length(gq_pts)
            
%             [i,j]
            
            zeta = gq_pts(i);
            weightz = gq_wts(i);            
            eta  = gq_pts(j);
            weighte = gq_wts(j);
                
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
                        
            [iis,jjs] = meshgrid(1:n);
%             [jjs,iis] = meshgrid(1:n);
            for k = 1:n_nodes
                q = elenode_order(k);
                
                % index for mesh
                ii = iis(q);
                jj = jjs(q);
                
                % zeta polynomial
                rz = zetas; rz(ii) = [];
                pz = poly(rz);
                pznorm = polyval(pz,zetas(ii));
                pz = pz/pznorm;
                
                % eta polynomial
                re = etas; re(jj) = [];
                pe = poly(re);
                penorm = polyval(pe,etas(jj));
                pe = pe/penorm;
                
                % calculate normalization coefficient
                pnorm = polyval(pz,zetas(ii))*polyval(pe,etas(jj));
                
                % shape function
                N(k) = polyval(pz,zeta)*polyval(pe,eta)/pnorm;
 
                % shape function derivatives
                dpzdz = pz(1:end-1).*[length(pz)-1:-1:1];
                dNdz(k) = polyval(dpzdz,zeta)*polyval(pe,eta)/pnorm;
                
                dpede = pe(1:end-1).*[length(pe)-1:-1:1];
                dNde(k) = polyval(pz,zeta)*polyval(dpede,eta)/pnorm;
            end
            
            J = [dNdz*xlocs,dNdz*ylocs;dNde*xlocs,dNde*ylocs];
            Rho = inv(J);
            
            dNdX = Rho(1,:)*[dNdz;dNde];
            dNdY = Rho(2,:)*[dNdz;dNde];
            
            jacob = det(J);
            jacobsum = jacobsum+jacob;
            
            % Shape Function Matrix
            Nmat = zeros(npn,ndof);
            
            for k = 1:n_nodes
                Nmat(:,(npn*(k-1)+1:npn*k)) = N(k)*eye(npn);
            end
                        
            % Element Mass Matrix
            m_add = weightz*weighte*jacob*(Nmat.'*rho_mat*Nmat);
            M_ele = M_ele+m_add;      
            
            % Element Bending Stiffness Matrix
            Bf = zeros(3,ndof);
            for k1 = 1:n_nodes
                Bf(:,(npn*(k1-1)+1:npn*k1)) = [0,dNdX(k1),0;...
                                               0,0,dNdY(k1);...
                                               0,dNdY(k1),dNdX(k1)];
            end
            
            k_add = weighte*weightz*jacob*(h^3/12)*(Bf'*Df*Bf);
            KB_ele = KB_ele + k_add;
            
            
            % Element Shear Stiffness Matrix
            Bc = zeros(2,ndof);
            for k1 = 1:n_nodes
                Bc(:,(npn*(k1-1)+1:npn*k1)) = [dNdX(k1), N(k1), 0;...
                                               dNdY(k1), 0,     N(k1)];
            end
            
            
            k_add = weighte*weightz*jacob*(h*alpha*Bc'*Dc*Bc);
            KS_ele_full = KS_ele_full+k_add;      
            
        end
    end   
    
    if jacobsum<0
        disp('negative jacobian')
    end
    %% Perform Reduced Integration of Shear Stiffness
    % Function lgwt generates gauss legendre points and weights
    [gq_pts,gq_wts] = lgwt(n-1,-1,1);
    
    for j = 1:length(gq_pts)
        for i = 1:length(gq_pts)

            zeta = gq_pts(i);
            weightz = gq_wts(i);            
            eta  = gq_pts(j);
            weighte = gq_wts(j);
                
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
                        
            [iis,jjs] = meshgrid(1:n);
%             [jjs,iis] = meshgrid(1:n);
            for k = 1:n_nodes
                q = elenode_order(k);
                
                % index for mesh
                ii = iis(q);
                jj = jjs(q);
                
                % zeta polynomial
                rz = zetas; rz(ii) = [];
                pz = poly(rz);
                pznorm = polyval(pz,zetas(ii));
                pz = pz/pznorm;
                
                % eta polynomial
                re = etas; re(jj) = [];
                pe = poly(re);
                penorm = polyval(pe,etas(jj));
                pe = pe/penorm;
                
                % calculate normalization coefficient
                pnorm = polyval(pz,zetas(ii))*polyval(pe,etas(jj));
                
                % shape function
                N(k) = polyval(pz,zeta)*polyval(pe,eta)/pnorm;
 
                % shape function derivatives
                dpzdz = pz(1:end-1).*[length(pz)-1:-1:1];
                dNdz(k) = polyval(dpzdz,zeta)*polyval(pe,eta)/pnorm;
                
                dpede = pe(1:end-1).*[length(pe)-1:-1:1];
                dNde(k) = polyval(pz,zeta)*polyval(dpede,eta)/pnorm;
            end
            
            J = [dNdz*xlocs,dNdz*ylocs;dNde*xlocs,dNde*ylocs];
            Rho = inv(J);
            
            dNdX = Rho(1,:)*[dNdz;dNde];
            dNdY = Rho(2,:)*[dNdz;dNde];
            
            jacob = det(J);
                        
            % Element Shear Stiffness Matrix
            Bc = zeros(2,ndof);
            for k1 = 1:n_nodes
                Bc(:,(npn*(k1-1)+1:npn*k1)) = [dNdX(k1), N(k1), 0;...
                                               dNdY(k1), 0,     N(k1)];
            end
            
            
            k_add = weighte*weightz*jacob*(h*alpha*Bc'*Dc*Bc);
            KS_ele_red = KS_ele_red+k_add;      
            
        end
    end    

% perturbation of reduced shear stiffness
% eps = 0.001;
KS_ele = KS_ele_red*(1-eps) + KS_ele_full*eps;

% Combined Element Stiffness
K_ele = KB_ele + KS_ele;

% K_ele
% pause

