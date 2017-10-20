function [B_ele] = element_BMatrix_3DE(n,kappa,...
    xlocs,ylocs,zlocs,elenode_order)

    n_kap = size(kappa,2);
    n_nodes = n^3;
    ndof = 3*n_nodes;
    
    [~,i_unsort_nodes] = sort(elenode_order);
    
    % Function lgwt generates gauss legendre points and weights
    node_pts = linspace(-1,1,n);

    
    % Preallocate arrays
    B_ele = zeros(6,ndof,n_nodes);

    % loop over gauss quadriture points to perform integration of shape 
    % functions, then create element stiffness and mass
    
    zetas = linspace(-1,1,n);
    etas = linspace(-1,1,n);
    ksis = linspace(-1,1,n);

    jacobsum = 0;
    for l = 1:n
        for j = 1:n
            for i = 1:n 
            
            
%                 
            
            zeta = node_pts(i);            
            eta  = node_pts(j);            
            ksi  = node_pts(l);
                            
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
            dNdk = zeros(1,n_nodes); 
            
            [jjs,iis,kks] = meshgrid([1:n],[1:n],[1:n]);
            for k = 1:n_nodes
                q = elenode_order(k);
%                 pause
                
                
                % index for mesh
                ii = iis(q);
                jj = jjs(q);
                kk = kks(q);
                                
                % zeta polynomial
                rz = zetas; rz(ii) = [];
                pz = poly(rz);
                
                % eta polynomial
                re = etas; re(jj) = [];
                pe = poly(re);
                
                % ksi polynomial
                rk = ksis; rk(kk) = [];
                pk = poly(rk);
                
                % calculate normalization coefficient
                pnorm = polyval(pz,zetas(ii))*...
                    polyval(pe,etas(jj))*...
                    polyval(pk,ksis(kk));
                
                % shape function
                N(k) = polyval(pz,zeta)*polyval(pe,eta)*polyval(pk,ksi)/pnorm;
 
                % shape function derivatives
                dpzdz = pz(1:end-1).*[length(pz)-1:-1:1];
                dNdz(k) = polyval(dpzdz,zeta)*...
                          polyval(pe,eta)*...
                          polyval(pk,ksi)/pnorm;
                
                dpede = pe(1:end-1).*[length(pe)-1:-1:1];
                dNde(k) = polyval(pz,zeta)*...
                    polyval(dpede,eta)*...
                    polyval(pk,ksi)/pnorm;
                
                dpkdk = pk(1:end-1).*[length(pk)-1:-1:1];
                dNdk(k) = polyval(pz,zeta)*...
                          polyval(pe,eta)*...
                          polyval(dpkdk,ksi)/pnorm;
            end
            
            J = [dNdz*xlocs,dNdz*ylocs,dNdz*zlocs;...
                dNde*xlocs,dNde*ylocs,dNde*zlocs;...
                dNdk*xlocs,dNdk*ylocs,dNdk*zlocs];
            
            jacob = det(J);
            jacobsum = jacobsum+jacob;
            Rho = inv(J);
            Rho_triple = [Rho,zeros(3,6);...
                    zeros(3),Rho,zeros(3);...
                    zeros(3,6),Rho];
            
%             % Shape Function Matrix
%             Nmat = zeros(3,ndof);
%             for k = 1:n_nodes
%                 Nmat(:,(3*(k-1)+1:3*k)) = N(k)*eye(3);
%             end
             
            % Element Stiffness Matrix
            B = zeros(6,ndof);
            
            % squish relates displacement derivatives to strains
            squish = zeros(9,6);
            squish([1,14,27,33,35,39,43,47,49]) = 1;
%             squish([1,14,27,29,31,42,43,48,52]) = 1;
            squish = squish';
            
            onelocs = zeros(9,3);
            onelocs([1,2,3,13,14,15,25,26,27]) = 1;
            
            
            for k1 = 1:n_nodes

                Bsmall = diag([dNdz(k1),dNde(k1),dNdk(k1),...
                                dNdz(k1),dNde(k1),dNdk(k1),...
                                dNdz(k1),dNde(k1),dNdk(k1)])*onelocs;
                
                    
                B(:,(3*k1-2:3*k1)) = squish*Rho_triple*Bsmall;
                
            end
                        
            node_count = (l-1)*n^2 + (j-1)*n + i;
            current_node = i_unsort_nodes(node_count);
            B_ele(:,:,current_node) = B;
            end

        end
    end















function [B_ele] = element_mass_stiffness_3DE(D,rho,n,coordinates)

n_nodes = n^3;
ndof = 3*n_nodes;

% Function lgwt generates gauss legendre points and weights
[gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);

% Preallocate arrays
B_ele = zeros(6,ndof,n_nodes);

% loop over gauss quadriture points to perform integration of shape 
% functions, then create element stiffness and mass

% compute coefficients of lagrange polynomials
L_p = zeros(n,n);
for i = 1:n        
    % polynomial roots
    root_x = linspace(-1,1,n);
    x_one = root_x(i);
    root_x(i) = [];        
    L_p(i,:) = poly(root_x);

    % normalize to unit value
    L_p(i,:) = L_p(i,:)/...
        polyval(L_p(i,:),x_one);        
end

% compute value of lagrange polynomials at GQ points
L_y = zeros(n,length(gq_pts));
for i = 1:n    
    L_y(i,:) = polyval(L_p(i,:),gq_pts);
end

% compute coefficients of lagrange polynomial derivative
dLdx_p = zeros(n,n-1);
for i = 1:n     
    dLdx_p(i,:) = L_p(i,1:end-1).*(n-1:-1:1);
end

% compute value of lagrange polynomials at GQ points
dLdx_y = zeros(n,length(gq_pts));
for i = 1:n    
    dLdx_y(i,:) = polyval(dLdx_p(i,:),gq_pts);
end

jacobsum = 0;


[iis,jjs,kks] = ndgrid(1:n,1:n,1:n);
for i = 1:length(gq_pts)
    for j = 1:length(gq_pts)
        for l = 1:length(gq_pts)

            % Gauss quadrature (GQ) weights for current GQ point
            weightz = gq_wts(i);  
            weighte = gq_wts(j);
            weightk = gq_wts(l);

            % preallocate shape function and shape function derivative
            % arrays
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
            dNdk = zeros(1,n_nodes); 

            for k = 1:n_nodes                

                % index for mesh
                ii = iis(k);
                jj = jjs(k);
                kk = kks(k);

                % shape function
                N(k) = L_y(ii,i)*L_y(jj,j)*L_y(kk,l);

                % shape function derivatives with respect to zeta, eta, and
                % ksi
                dNdz(k) = dLdx_y(ii,i)*L_y(jj,j)*L_y(kk,l);
                dNde(k) = L_y(ii,i)*dLdx_y(jj,j)*L_y(kk,l);
                dNdk(k) = L_y(ii,i)*L_y(jj,j)*dLdx_y(kk,l);
            end

            % Jacobian matrix            
            J = [dNdz;dNde;dNdk]*coordinates;

            % Jacobian determinant
            jacob = det(J);
            jacobsum = jacobsum+jacob;


            dNdX = J\[dNdz;dNde;dNdk];

            % Shape Function Matrix
            Nmat = zeros(3,ndof);
            for k = 1:n_nodes
                Nmat(:,(3*(k-1)+1:3*k)) = N(k)*eye(3);
            end

            % Element Mass Matrix
            m_add = weightz*weighte*weightk*jacob*rho*(Nmat.'*Nmat);
            M_ele(:,:) = M_ele(:,:)+m_add;      

            % Element Stiffness Matrix
            B = zeros(6,ndof);

            B(1,1:3:end) = dNdX(1,:);
            B(2,2:3:end) = dNdX(2,:);
            B(3,3:3:end) = dNdX(3,:);
            B(4,2:3:end) = dNdX(3,:);
            B(4,3:3:end) = dNdX(2,:);
            B(5,1:3:end) = dNdX(3,:);
            B(5,3:3:end) = dNdX(1,:);
            B(6,1:3:end) = dNdX(2,:);
            B(6,2:3:end) = dNdX(1,:);

            k_add = weightz*weighte*weightk*jacob*(B'*D*B);
            K_ele = K_ele+k_add;
        end
    end
end
    
if jacobsum<0
    disp('negative jacobian')
end