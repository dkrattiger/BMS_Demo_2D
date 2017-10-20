function [H_ele,S_ele] = ...
    element_matrices_3D_Electronic(n,coordinates,V)

n_nodes = n^3;
ndof = n_nodes;

% Function lgwt generates gauss legendre points and weights
[gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);

% Preallocate arrays
H_ele = zeros(ndof,ndof);
S_ele = zeros(ndof,ndof);

% form Lagrange polynomial templates
lgnd = zeros(n,n);
for i = 1:n
    r = linspace(-1,1,n);
    rn = r(i);
    r(i) = [];
    
    % create a Lagrange polynomial that 
    lgnd(i,:) = poly(r)/polyval(poly(r),rn);
end

% evaluate Lagrange polynomial templates at node/gq pt combinations
lgnd_eval = zeros(n,length(gq_pts));
for i = 1:n
    lgnd_eval(i,:) = polyval(lgnd(i,:),gq_pts);
end

% form Lagrange polynomial derivative templates
d_lgnd = zeros(n,n-1);
for i = 1:n
    d_lgnd(i,:) = lgnd(i,1:end-1).*(n-1:-1:1);
end

% evaluate Lagrange polynomial templates at node/gq pt combinations
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

%             [iis,jjs,kks] = meshgrid([1:n],[1:n],[1:n]);
            [iis,jjs,kks] = ndgrid([1:n],[1:n],[1:n]);
            
            % loop through nodes and compute shape functions and shape
            % function derivatives for current gq point
            for k = 1:n_nodes

                % index for mesh
                ii = iis(k);
                jj = jjs(k);
                kk = kks(k);

                % shape function (is the product of 3 legendre polynomials
                % in different coordinates
                N(k) = lgnd_eval(ii,i)*lgnd_eval(jj,j)*lgnd_eval(kk,l);
                
                % shape function derivatives
                dNdz(k) = d_lgnd_eval(ii,i) *   lgnd_eval(jj,j)     *   lgnd_eval(kk,l);
                dNde(k) = lgnd_eval(ii,i)   *   d_lgnd_eval(jj,j)   *   lgnd_eval(kk,l);
                dNdk(k) = lgnd_eval(ii,i)   *   lgnd_eval(jj,j)     *   d_lgnd_eval(kk,l);
            end

%             J = [dNdz*xlocs,dNdz*ylocs,dNdz*zlocs;...
%                  dNde*xlocs,dNde*ylocs,dNde*zlocs;...
%                  dNdk*xlocs,dNdk*ylocs,dNdk*zlocs];
             
            J = [dNdz;dNde;dNdk]*coordinates;

            jacob = det(J);
            jacobsum = jacobsum+jacob;
%             Rho = inv(J);

            % Element Mass Matrix
            S_add = weightz*weighte*weightk*jacob*(N.'*N);  
            S_ele(:,:) = S_ele(:,:)+S_add;
            
            % strain-displacement matrix
            B = J\[dNdz;dNde;dNdk];

            % potential
            % =========         
            if isa(V,'function_handle')

                % determine realspace position vector to compute
                % pseudopotential
%                 rvec = zeros(3,1);
%                 rvec(1,1) = N*xlocs;
%                 rvec(2,1) = N*ylocs;
%                 rvec(3,1) = N*zlocs;
                rvec = (N*coordinates)';

                Vloc = V(rvec);
            else
                Vloc = V;   
            end

            hbar = 1.05457e-34;
            me = 9.10938e-31;
            hbar2_me = hbar^2/me;

            hbar2_me = 1;


            H_add_V = weightz*weighte*weightk*jacob*Vloc*(N'*N);
            H_add_0 = weightz*weighte*weightk*jacob*(B'*B);

            H_ele  = H_ele  + hbar2_me*H_add_0  + H_add_V;
        end
    end
end

if jacobsum<0
    disp('negative jacobian')
end