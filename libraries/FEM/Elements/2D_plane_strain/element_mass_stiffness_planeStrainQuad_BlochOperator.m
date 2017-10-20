function [K0_ele,Kx_ele,Ky_ele,Kxx_ele,Kyy_ele,Kxy_ele,M_ele] ...
    = element_mass_stiffness_planeStrainQuad_BlochOperator(D,rho,n,kappa,...
    xlocs,ylocs,elenode_order)

    n_kap = size(kappa,2);
    n_nodes = n^2;
    ndof = 2*n_nodes;
    
    % Function lgwt generates gauss legendre points and weights
    [gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);

    
    % Preallocate arrays
    K0_ele = zeros(ndof);
    Kx_ele = zeros(ndof);
    Ky_ele = zeros(ndof);
    Kxx_ele = zeros(ndof);
    Kyy_ele = zeros(ndof);
    Kxy_ele = zeros(ndof);
    K_ele = zeros(ndof,ndof,n_kap);
        M_ele = zeros(ndof);

    % loop over gauss quadriture points to perform integration of shape 
    % functions, then create element stiffness and mass
    
    zetas = linspace(-1,1,n);
    etas = linspace(-1,1,n);

    jacobsum = 0;
    for j = 1:length(gq_pts)
        for i = 1:length(gq_pts)
%         i = ij(1,q);
%         j = ij(2,q);

            zeta = gq_pts(i);
            weightz = gq_wts(i);            
            eta  = gq_pts(j);
            weighte = gq_wts(j);
                
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
                        
            [jjs,iis] = meshgrid(1:n);
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

            jacob = det(J);
            jacobsum = jacobsum+jacob;
            Rho_double = [inv(J),zeros(2);zeros(2),inv(J)];
            
            % Shape Function Matrix
            Nmat = zeros(2,ndof);
            for k = 1:n_nodes
                Nmat(:,(2*k-1:2*k)) = N(k)*eye(2);
            end
            
            % Element Mass Matrix
            m_add = weightz*weighte*jacob*rho*(Nmat.'*Nmat);
            M_ele(:,:) = M_ele(:,:)+m_add;      
            
            % Element Stiffness Matrix
            B = zeros(3,ndof);
            squish = [1,0,0,0;0,0,0,1;0,1,1,0];
            for k1 = 1:n_nodes
                
                Bsmall = [dNdz(k1),0;dNde(k1),0;0,dNdz(k1);0,dNde(k1)];
                B(:,(2*k1-1:2*k1)) = squish*Rho_double*Bsmall;
            end
            
%             Bk = zeros(3,ndof);
%             for k2 = 1:n_kap
%                 kx = kappa(1,k2);
%                 ky = kappa(2,k2);
%                 
%                 for k3 = 1:n_nodes
%                     Bksmall = [kx*N(k3),0;ky*N(k3),0;0,kx*N(k3);0,ky*N(k3)];
%                     Bk(:,(2*k3-1:2*k3)) = 1i*squish*Bksmall;
%                 end
% 
% %                 size(1i*squish*Bksmall)
% 
%                 Btilde = B+Bk;
%                 k_add = weighte*weightz*jacob*(Btilde'*D*Btilde);
%                 K_ele(:,:,k2) = K_ele(:,:,k2)+k_add;
%             end
            Bkx = zeros(3,ndof);
            Bky = zeros(3,ndof);
            for k3 = 1:n_nodes
                kx = 1;
                ky = 1;
                Bkxsmall = [kx*N(k3),0;0,0;0,kx*N(k3);0,0];
                Bkx(:,(2*k3-1:2*k3)) = 1i*squish*Bkxsmall;
                Bkysmall = [0,0;ky*N(k3),0;0,0;0,ky*N(k3)];
                Bky(:,(2*k3-1:2*k3)) = 1i*squish*Bkysmall;
            end
            
            k0_add = weighte*weightz*jacob*(B'*D*B);
            K0_ele = K0_ele+k0_add;

            kxx_add = weighte*weightz*jacob*(Bkx'*D*Bkx);
            Kxx_ele = Kxx_ele+kxx_add;
            
            kyy_add = weighte*weightz*jacob*(Bky'*D*Bky);
            Kyy_ele = Kyy_ele+kyy_add;
            
            kx_add = weighte*weightz*jacob*(B'*D*Bkx+Bkx'*D*B);
            Kx_ele = Kx_ele+kx_add;
            
            ky_add = weighte*weightz*jacob*(B'*D*Bky+Bky'*D*B);
            Ky_ele = Ky_ele+ky_add;
            
            kxy_add = weighte*weightz*jacob*(Bkx'*D*Bky+Bky'*D*Bkx);
            Kxy_ele = Kxy_ele+kxy_add;

        end
    end    
    
if jacobsum<0
    disp('negative jacobian')
end