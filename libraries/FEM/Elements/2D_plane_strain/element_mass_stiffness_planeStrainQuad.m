function [K_ele,M_ele] ...
    = element_mass_stiffness_planeStrainQuad(D,rho,n,xlocs,ylocs)
    
% Function "element_mass_stiffness_planeStrainQuad" generates the
% element mass and stiffness matrices for an (n x n)-node plane strain
% finite element. 
% 
% Inputs:
% D         =   element elasticity matrix
% rho       =   element density
% n         =   number of nodes along element edge (element order + 1)
% xlocs,ylocs     
%           =   x, y coordinates of element nodes 
%               (note: pay attention to the ordering)

n_nodes = n^2;
ndof = 2*n_nodes;

% Function lgwt generates gauss legendre points and weights
[gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);

%% Form and evaluate 1D shape functions
% ======================================================================= %
% form 1D lagrange polynomials 
% (these get tensor multiplied to get to 2D)
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


% Preallocate arrays
K_ele = zeros(ndof);
M_ele = zeros(ndof);


%% Loop Over Gauss Points
% ======================================================================= %
% loop over gauss quadriture points to perform integration of shape 
% functions, then create element stiffness and mass

jacobsum = 0; % jacobian sum
for j = 1:length(gq_pts)
    for i = 1:length(gq_pts)

        % Gauss quadrature weights
        weightz = gq_wts(i);
        weighte = gq_wts(j);

        % preallocate shape function and shape function derivative
        % vectors
        N = zeros(1,n_nodes);
        dNdz = zeros(1,n_nodes);
        dNde = zeros(1,n_nodes);

        % I can't remember why I used "meshgrid" because I usually
        % prefer "ndgrid"
        [jjs,iis] = meshgrid(1:n);
        %[jjs,iis] = ndgrid(1:n);

        for k = 1:n_nodes
            q = k;

            % index for mesh
            ii = iis(q);
            jj = jjs(q);

            % shape function
            N(k) = lgnd_eval(ii,i)*lgnd_eval(jj,j);

            % shape function derivatives
            dNdz(k) = d_lgnd_eval(ii,i) *   lgnd_eval(jj,j);
            dNde(k) = lgnd_eval(ii,i)   *   d_lgnd_eval(jj,j);
        end      

        % Jacobian matrix
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
        k0_add = weighte*weightz*jacob*(B'*D*B);
        K_ele = K_ele+k0_add;
    end
end    
    
if jacobsum<0
    disp('negative jacobian')
end