function [B_ele] = element_BMatrix_3DE(n,coordinates)

n_nodes = n^3;
ndof = 3*n_nodes;

% Generate isoparametric node locations
node_pts = linspace(-1,1,n);

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

% compute value of lagrange polynomials at each node (this should be
% perfunctory due to the definition of lagrange polynomials)
L_y = zeros(n,length(node_pts));
for i = 1:n    
    L_y(i,:) = polyval(L_p(i,:),node_pts);
end

% compute coefficients of lagrange polynomial derivatives
dLdx_p = zeros(n,n-1);
for i = 1:n     
    dLdx_p(i,:) = L_p(i,1:end-1).*(n-1:-1:1);
end

% compute value of lagrange polynomials at each node
dLdx_y = zeros(n,length(node_pts));
for i = 1:n    
    dLdx_y(i,:) = polyval(dLdx_p(i,:),node_pts);
end


% Preallocate arrays
B_ele = zeros(6,ndof,n_nodes);

% loop over gauss quadriture points to perform integration of shape 
% functions, then create element stiffness and mass
[iis,jjs,kks] = ndgrid([1:n],[1:n],[1:n]);
for k = 1:n
    for j = 1:n
        for i = 1:n 

            % preallocate shape function and shape function derivative
            % arrays
            N = zeros(1,n_nodes);
            dNdz = zeros(1,n_nodes);
            dNde = zeros(1,n_nodes);
            dNdk = zeros(1,n_nodes); 

            for q = 1:n_nodes                

                % index for mesh
                ii = iis(q);
                jj = jjs(q);
                kk = kks(q);

                % shape function derivatives with respect to zeta, eta, and
                % ksi
                dNdz(q) = dLdx_y(ii,i)*L_y(jj,j)*L_y(kk,k);
                dNde(q) = L_y(ii,i)*dLdx_y(jj,j)*L_y(kk,k);
                dNdk(q) = L_y(ii,i)*L_y(jj,j)*dLdx_y(kk,k);
            end

            % Jacobian matrix            
            J = [dNdz;dNde;dNdk]*coordinates;

            % realspace shapefunction derivatives
            dNdX = J\[dNdz;dNde;dNdk];
            node_count = (k-1)*n^2 + (j-1)*n + i;

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

            B_ele(:,:,node_count) = B;
        end
    end
end