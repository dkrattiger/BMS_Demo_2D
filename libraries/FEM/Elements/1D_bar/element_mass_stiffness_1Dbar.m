function [K_ele,M_ele] = element_mass_stiffness(E,A,rho,n,kappa,xlocs)

    n_kap = length(kappa);

    % Function lgwt generates gauss legendre points and weights
    [gq_pts,gq_wts] = lgwt(n,-1,1);


    % Preallocate arrays
    K_ele = zeros(n,n,n_kap);
    M_ele = zeros(n,n);
    N = zeros(1,n);
    dNdz = zeros(1,n);

    % loop over gauss quadriture points to perform integration of shape 
    % functions, then create element stiffness and mass
    for j = 1:length(gq_pts)
        zeta = gq_pts(j);
        weight = gq_wts(j);
        zetas = linspace(-1,1,n);

        % loop over shape functions
        for k = 1:n;

            % polynomial roots
            r = zetas;
            r(k) = [];

            % create polynomial
            p = poly(r);

            % normalize polynomial
            mult = 1/polyval(p,zetas(k));
            p = p*mult;
            dpdz = p(1:end-1).*[length(p)-1:-1:1];

            % evaluate polynomials at gq point
            N(k) = polyval(p,zeta);
            dNdz(k) =  polyval(dpdz,zeta);

        end
        jacob = dNdz*xlocs';

        dNdx = dNdz/jacob;

        m_add = weight*jacob*rho*A*(N.'*N);
        M_ele(:,:) = M_ele(:,:)+m_add;

        B = dNdx;
        for k = 1:n_kap
            Bk = 1i*kappa(k)*N;

            Btilde =  B+Bk;
            k_add = weight*jacob*E*A*(Btilde'*Btilde);
            K_ele(:,:,k) = K_ele(:,:,k)+k_add;
        end
    end    