function [K_ele,M_ele] ...
    = element_mass_stiffness_mindlin_bloch_op(Df,Dc,rho,h,n,...
    xlocs,ylocs,elenode_order,eps)
    
%     elenode_order = 1:n^2;

    % different kxky combinations
    field_names = {'s0','sx','sy','sxx','sxy','syy'};
    
    % shear parameter
    alpha = 5/6; 
    
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
%     KB_ele = zeros(ndof);
%     KS_ele_full = zeros(ndof);
%     KS_ele_red = zeros(ndof);
%     M_ele = zeros(ndof);
    for k = 1:length(field_names)
        KB_ele.(field_names{k}) = zeros(ndof);  
        KS_ele_full.(field_names{k}) = zeros(ndof);  
        KS_ele_red.(field_names{k}) = zeros(ndof);  
        
    end
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
                        
%             [jjs,iis] = meshgrid(1:n);
            [iis,jjs] = meshgrid(1:n);
            for k = 1:n_nodes
                q = elenode_order(k);
%                 q=k;
                
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
            Bf0 = zeros(3,ndof);
            Bfkx = zeros(3,ndof);
            Bfky = zeros(3,ndof);
            for k1 = 1:n_nodes
                Bf0(:,(npn*(k1-1)+1:npn*k1)) = [0,  dNdX(k1),   0;...
                                                0,  0,          dNdY(k1);...
                                                0,  dNdY(k1),   dNdX(k1)];
                                            
                Bfkx(:,(npn*(k1-1)+1:npn*k1)) = [0, 1i*N(k1),   0;...
                                                 0, 0,          0;...
                                                 0, 0,          1i*N(k1)];
                                            
                Bfky(:,(npn*(k1-1)+1:npn*k1)) = [0, 0,          0;...
                                                 0, 0,          1i*N(k1);...
                                                 0, 1i*N(k1),   0];
            end
            
%             k_add = weighte*weightz*jacob*(h^3/12)*(Bf0'*Df*Bf0);
            k_add.s0 = weighte*weightz*jacob*(h^3/12)*(Bf0'*Df*Bf0);
            k_add.sx = weighte*weightz*jacob*(h^3/12)*(Bfkx'*Df*Bf0 + Bf0'*Df*Bfkx);
            k_add.sy = weighte*weightz*jacob*(h^3/12)*(Bfky'*Df*Bf0 + Bf0'*Df*Bfky);
            k_add.sxx = weighte*weightz*jacob*(h^3/12)*(Bfkx'*Df*Bfkx);
            k_add.sxy = weighte*weightz*jacob*(h^3/12)*(Bfkx'*Df*Bfky + Bfky'*Df*Bfkx);
            k_add.syy = weighte*weightz*jacob*(h^3/12)*(Bfky'*Df*Bfky);
%             KB_ele = KB_ele + k_add;
            for k = 1:length(field_names)
                KB_ele.(field_names{k}) = ...
                    KB_ele.(field_names{k})+k_add.(field_names{k});      
            end
            
%             % Element Shear Stiffness Matrix
%             Bc0 = zeros(2,ndof);
%             for k1 = 1:n_nodes
%                 Bc0(:,(npn*(k1-1)+1:npn*k1)) = [dNdX(k1), -N(k1), 0;...
%                                                dNdY(k1), 0,     -N(k1)];
%             end
%             
%             
%             k_add = weighte*weightz*jacob*(h*alpha*Bc0'*Dc*Bc0);
%             KS_ele_full = KS_ele_full+k_add;     

            % Element Shear Stiffness Matrix
            Bc0 = zeros(2,ndof);
            Bckx = zeros(2,ndof);
            Bcky = zeros(2,ndof);
            for k1 = 1:n_nodes
                Bc0(:,(npn*(k1-1)+1:npn*k1)) = [dNdX(k1),   N(k1), 0;...
                                                dNdY(k1),   0,     N(k1)];
                
                Bckx(:,(npn*(k1-1)+1:npn*k1)) = [N(k1)*1i,  0,      0;...
                                                 0,         0,      0];
                Bcky(:,(npn*(k1-1)+1:npn*k1)) = [0,         0,      0;...
                                                 N(k1)*1i,  0,      0];
            end
            
            k_add.s0 = weighte*weightz*jacob*h*alpha*(Bc0'*Dc*Bc0);
            k_add.sx = weighte*weightz*jacob*h*alpha*(Bckx'*Dc*Bc0 + Bc0'*Dc*Bckx);
            k_add.sy = weighte*weightz*jacob*h*alpha*(Bcky'*Dc*Bc0 + Bc0'*Dc*Bcky);
            k_add.sxx = weighte*weightz*jacob*h*alpha*(Bckx'*Dc*Bckx);
            k_add.sxy = weighte*weightz*jacob*h*alpha*(Bckx'*Dc*Bcky + Bcky'*Dc*Bckx);
            k_add.syy = weighte*weightz*jacob*h*alpha*(Bcky'*Dc*Bcky);
            
            for k = 1:length(field_names)
                KS_ele_full.(field_names{k}) = ...
                    KS_ele_full.(field_names{k})+k_add.(field_names{k});      
            end
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
                        
%             [jjs,iis] = meshgrid(1:n);
            [iis,jjs] = meshgrid(1:n);
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
            Bc0 = zeros(2,ndof);
            Bckx = zeros(2,ndof);
            Bcky = zeros(2,ndof);
            for k1 = 1:n_nodes
                Bc0(:,(npn*(k1-1)+1:npn*k1)) = [dNdX(k1), N(k1), 0;...
                                               dNdY(k1), 0,     N(k1)];
                
                Bckx(:,(npn*(k1-1)+1:npn*k1)) = [N(k1)*1i, 0, 0;...
                                                 0, 0, 0];
                Bcky(:,(npn*(k1-1)+1:npn*k1)) = [0, 0, 0;...
                                                 N(k1)*1i, 0, 0];
            end
            
            
            k_add.s0 = weighte*weightz*jacob*h*alpha*(Bc0'*Dc*Bc0);
            k_add.sx = weighte*weightz*jacob*h*alpha*(Bckx'*Dc*Bc0 + Bc0'*Dc*Bckx);
            k_add.sy = weighte*weightz*jacob*h*alpha*(Bcky'*Dc*Bc0 + Bc0'*Dc*Bcky);
            k_add.sxx = weighte*weightz*jacob*h*alpha*(Bckx'*Dc*Bckx);
            k_add.sxy = weighte*weightz*jacob*h*alpha*(Bckx'*Dc*Bcky + Bcky'*Dc*Bckx);
            k_add.syy = weighte*weightz*jacob*h*alpha*(Bcky'*Dc*Bcky);
            
            for k = 1:length(field_names)
                KS_ele_red.(field_names{k}) = ...
                    KS_ele_red.(field_names{k})+k_add.(field_names{k});      
            end
            
        end
    end    


% eps = 0.001;
for k = 1:length(field_names)
    % perturbation of reduced shear stiffness
    KS_ele.(field_names{k}) = ...
            KS_ele_red.(field_names{k})*(1-eps)+...
            KS_ele_full.(field_names{k})*eps;
        
    % Combined Element Stiffness
    K_ele.(field_names{k}) = KB_ele.(field_names{k}) + KS_ele.(field_names{k});    
end



