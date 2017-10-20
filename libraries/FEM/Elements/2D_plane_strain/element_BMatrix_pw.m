function B_ele = element_BMatrix_pw(n,kappa,xlocs,ylocs,elenode_order)

n_kap = size(kappa,2);
n_nodes = n^2;
ndof = 2*n_nodes;
[~,i_unsort_nodes] = sort(elenode_order);

% Function lgwt generates gauss legendre points and weights
[gq_pts,gq_wts] = lgwt((n-1)*2,-1,1);


% Preallocate arrays
B_ele = zeros(3,ndof,n_nodes);

% loop over gauss quadriture points to perform integration of shape 
% functions, then create element stiffness and mass

zetas = linspace(-1,1,n);
etas = linspace(-1,1,n);

jacobsum = 0;
for j = 1:length(gq_pts)
    for i = 1:length(gq_pts)

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

%         % Shape Function Matrix
%         Nmat = zeros(2,ndof);
%         for k = 1:n_nodes
%             Nmat(:,(2*k-1:2*k)) = N(k)*eye(2);
%         end

%         % Element Mass Matrix
%         m_add = weightz*weighte*jacob*rho*(Nmat.'*Nmat);
%         M_ele(:,:) = M_ele(:,:)+m_add;      

        % Element Stiffness Matrix
        B = zeros(3,ndof);
        squish = [1,0,0,0;0,0,0,1;0,1,1,0];
        for k1 = 1:n_nodes

            Bsmall = [dNdz(k1),0;dNde(k1),0;0,dNdz(k1);0,dNde(k1)];
            B(:,(2*k1-1:2*k1)) = squish*Rho_double*Bsmall;
        end

                        
        node_count = (j-1)*n + i;
        
        current_node = i_unsort_nodes(node_count);
        B_ele(:,:,current_node) = B;
    end
end    
    
if jacobsum<0
    disp('negative jacobian')
end