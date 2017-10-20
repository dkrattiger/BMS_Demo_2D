function [K_ele,M_ele] = element_mass_stiffness_EB1D(E,rho,I,A,kappa,...
        xlocs,options)

% Default options
defaults.elementType = 'standard'; 
    
% fill in unspecified options fields with the default values
options = setstructfields(defaults,options);

n_kap = size(kappa,2);
l = abs(xlocs(1)-xlocs(2));

K_ele = zeros(4,4,n_kap); 
if strcmpi(options.elementType,'unconjugatedBlochOperator')
    M_ele = zeros(4,4,n_kap); 
else
    M_ele = zeros(4,4); 
end

for j = 1:n_kap
    k = kappa(j);

    % polynomial representations of the shape functions
    pH = [2/(l^3),    -3/(l^2),  0,      1;...
          1/(l^2),    -2/(l),    1,      0;...
          -2/(l^3),    3/(l^2),  0,      0;...
          1/(l^2),    -1/(l),   0,      0];
    pH_x = pH(:,1:end-1).*(ones(4,1)*((size(pH,2)-1):-1:1));
    pH_xx = pH_x(:,1:end-1).*(ones(4,1)*((size(pH_x,2)-1):-1:1));

    % gauss-quadrature points
    nqp = 5;
    [qx,qw] = lgwt(nqp,0,l);

    n_shapefun = size(pH,1);
    pHeval = zeros(n_shapefun,nqp);
    pH_xeval = zeros(n_shapefun,nqp);
    pH_xxeval = zeros(n_shapefun,nqp);
    for i_shapefun = 1:n_shapefun
        pHeval(i_shapefun,:) = polyval(pH(i_shapefun,:),qx);
        pH_xeval(i_shapefun,:) = polyval(pH_x(i_shapefun,:),qx);   
        pH_xxeval(i_shapefun,:) = polyval(pH_xx(i_shapefun,:),qx);              
    end

    for i_qp = 1:nqp

        % gauss quadrature point
        xqp = qx(i_qp);
        wqp = qw(i_qp);

        % when integrating the exponential component of the
        % shape functions, it is important to use the absolute
        % (not element) value of x
        x = xlocs(1)+xqp;

        % traditional EB shape functions and their derivatives evaluated at
        % the current qp point
        H = pHeval(:,i_qp).';
        H_x= pH_xeval(:,i_qp).';
        H_xx= pH_xxeval(:,i_qp).';

        % exponential terms do not cancel when evaluating unconjugated
        % version
        if strcmpi(options.elementType,'unconjugatedBlochOperator')
            N = H*exp(1i*k*x);
            B = (H_xx + 2*1i*k*H_x - k^2*H)*exp(1i*k*x);
            K_ele(:,:,j) = K_ele(:,:,j) + (wqp*E*I) * (B.'*B);
            M_ele(:,:,j) = M_ele(:,:,j) + (wqp*rho*A) * (N.'*N);
        else
            N = H;
            B = (H_xx + 2*1i*k*H_x - k^2*H);   
            K_ele(:,:,j) = K_ele(:,:,j) + (wqp*E*I) * (B'*B);
            if j==1
                M_ele = M_ele + (wqp*rho*A) * (N'*N);
            end
        end


    end
end
