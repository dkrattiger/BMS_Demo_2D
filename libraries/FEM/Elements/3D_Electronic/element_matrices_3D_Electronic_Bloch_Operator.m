function [H_ele_0,H_ele_kx,H_ele_ky,H_ele_kz,S_ele] = element_matrices_3D_Electronic_Bloch_Operator(n,kappa,...
    xlocs,ylocs,zlocs,V,elenode_order)


tic
n_kap = size(kappa,2);
n_nodes = n^3;
ndof = n_nodes;

% Function lgwt generates gauss legendre points and weights
[gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);

% Preallocate arrays
H_ele_0 = zeros(ndof,ndof);
H_ele_kx = zeros(ndof,ndof);
H_ele_kz = zeros(ndof,ndof);
H_ele_ky = zeros(ndof,ndof);
S_ele = zeros(ndof,ndof);

% form legendre polynomial templates
lgnd = zeros(n,n);
for i = 1:n
    r = linspace(-1,1,n);
    rn = r(i);
    r(i) = [];
    
    % create a legendre polynomial that 
    lgnd(i,:) = poly(r)/polyval(poly(r),rn);
end

% evaluate legendre polynomial templates at node/gq pt combinations
lgnd_eval = zeros(n,length(gq_pts));
for i = 1:n
    lgnd_eval(i,:) = polyval(lgnd(i,:),gq_pts);
end

% form legendre polynomial derivative templates
d_lgnd = zeros(n,n-1);
for i = 1:n
    d_lgnd(i,:) = lgnd(i,1:end-1).*(n-1:-1:1);
end

% evaluate legendre polynomial templates at node/gq pt combinations
d_lgnd_eval = zeros(n,length(gq_pts));
for i = 1:n
    d_lgnd_eval(i,:) = polyval(d_lgnd(i,:),gq_pts);
end

% loop over gauss quadriture points to perform integration of shape 
% functions, then create element stiffness and mass
jacobsum = 0;
for i = 1:length(gq_pts) 
    for j = 1:length(gq_pts)
        for l = 1:length(gq_pts)

            % Gauss-legendre weights in each direction
            weightz = gq_wts(i);  
            weighte = gq_wts(j);
            weightk = gq_wts(l);

            % preallocate shape functions
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
            dNdk = zeros(1,n_nodes); 

            [iis,jjs,kks] = meshgrid([1:n],[1:n],[1:n]);
            
            % loop through nodes and compute shape functions and shape
            % function derivatives for current gq point
            for k = 1:n_nodes
                q = elenode_order(k);

                % index for mesh
                ii = iis(q);
                jj = jjs(q);
                kk = kks(q);

                % shape function (is the product of 3 legendre polynomials
                % in different coordinates
                N(k) = lgnd_eval(ii,i)*lgnd_eval(jj,j)*lgnd_eval(kk,l);
                
                % shape function derivatives
                dNdz(k) = d_lgnd_eval(ii,i) *   lgnd_eval(jj,j)     *   lgnd_eval(kk,l);
                dNde(k) = lgnd_eval(ii,i)   *   d_lgnd_eval(jj,j)   *   lgnd_eval(kk,l);
                dNdk(k) = lgnd_eval(ii,i)   *   lgnd_eval(jj,j)     *   d_lgnd_eval(kk,l);
            end

            J = [dNdz*xlocs,dNdz*ylocs,dNdz*zlocs;...
                 dNde*xlocs,dNde*ylocs,dNde*zlocs;...
                 dNdk*xlocs,dNdk*ylocs,dNdk*zlocs];

            jacob = det(J);
            jacobsum = jacobsum+jacob;
%             Rho = inv(J);

            % Element Mass Matrix
            S_add = weightz*weighte*weightk*jacob*(N.'*N);  
            S_ele(:,:) = S_ele(:,:)+S_add;      

            % Element Stiffness Matrix
%             B = zeros(6,ndof);

            % squish relates displacement derivatives to strains
%             squish = zeros(9,6);
%             squish([1,14,27,33,35,39,43,47,49]) = 1;
%             squish([1,14,27,29,31,42,43,48,52]) = 1;
%             squish = squish';

%             onelocs = zeros(9,3);
%             onelocs([1,2,3,13,14,15,25,26,27]) = 1;

%             B = zeros(3,ndof);
%             for k1 = 1:n_nodes
%                 B(:,k1) = Rho*[dNdz(k1);dNde(k1);dNdk(k1)];
%             end
            
            B = J\[dNdz;dNde;dNdk];
%             norm(B-B2)
%             pause

            % potential
            % =========
            if isa(V,'function_handle')

                % determine realspace position vector to compute
                % pseudopotential
                rvec = zeros(3,1);
                rvec(1,1) = N*xlocs;
                rvec(2,1) = N*ylocs;
                rvec(3,1) = N*zlocs;

                Vloc = V(rvec);
            else
                Vloc = V;
            end

                % compute pseudopotential
%                     V = Si_Pseudo(rvec,a0);

%                 for k2 = 1:n_kap
%                     kx = kappa(1,k2);
%                     ky = kappa(2,k2);
%                     kz = kappa(3,k2);
%                     
%                     ksquared = kx^2+ky^2+kz^2;

%                     H_add = weightz*weighte*weightk*jacob*(B'*B - ...
%                         2*sqrt(-1)*N'*kappa(:,k2)'*B + (V+ksquared)*N'*N);
%                     H_ele(:,:,k2) = H_ele(:,:,k2)+ H_add;

                hbar = 1.05457e-34;
                me = 9.10938e-31;
                hbar2_me = hbar^2/me;

                hbar2_me = 1;


                H_add_V = weightz*weighte*weightk*jacob*Vloc*(N'*N);
                H_add_0 = weightz*weighte*weightk*jacob*(B'*B);
                H_add_kx = -weightz*weighte*weightk*jacob*...
                    2i*(N'*[1,0,0]*B);
                H_add_ky = -weightz*weighte*weightk*jacob*...
                    2i*(N'*[0,1,0]*B);
                H_add_kz = -weightz*weighte*weightk*jacob*...
                    2i*(N'*[0,0,1]*B);

                H_ele_0  = H_ele_0  + hbar2_me*H_add_0   + H_add_V;
                H_ele_kx = H_ele_kx + hbar2_me*H_add_kx;
                H_ele_ky = H_ele_ky + hbar2_me*H_add_ky;
                H_ele_kz = H_ele_kz + hbar2_me*H_add_kz;

%                     H_ele_0  = H_ele_0  + hbar2_me*H_add_0;
%                     H_ele_kx = H_ele_kx + sqrt(hbar2_me)*H_add_kx;
%                     H_ele_ky = H_ele_ky + sqrt(hbar2_me)*H_add_ky;
%                     H_ele_kz = H_ele_kz + sqrt(hbar2_me)*H_add_kz;
%                 end
        end
    end
end
toc
if jacobsum<0
        disp('negative jacobian')
    end