function [K,M] = master_mass_stiffness_EB1D(kappa,Model_Info,options);


approach = 'numeric';
switch approach
    case 'numeric'
        % Default options
        defaults.elementType = 'standard'; 

        % fill in unspecified options fields with the default values
        options = setstructfields(defaults,options);

        % Mesh info
        Es              = Model_Info.Es;
        rhos            = Model_Info.rhos;
        Is              = Model_Info.Is;
        As              = Model_Info.As;
        pattern         = Model_Info.pattern;
        n_dof           = Model_Info.n_dof;
        elenodes        = Model_Info.elenodes;
        xlocs           = Model_Info.xlocs;

        Lx = max(xlocs)-min(xlocs);


        n_kap  = size(kappa,2);
        n_eles = size(elenodes,1);

        K = zeros(n_dof,n_dof,n_kap);
        if strcmpi(options.elementType,'unconjugatedBlochOperator')
            M = zeros(n_dof,n_dof,n_kap);
        else
            M = zeros(n_dof,n_dof);
        end

        for i = 1:n_eles
            i
            dofs_in_ele = zeros(4,1);
            for j = 1:length(elenodes(i,:))
                en = elenodes(i,j);
                dofs_in_ele(2*(j-1)+1:2*j) = [2*(en-1)+1:2*en];
            end

            rho = rhos{pattern(i)};
            E   = Es{pattern(i)};
            I   = Is{pattern(i)};
            A   = As{pattern(i)};

            [K_ele,M_ele] = element_mass_stiffness_EB1D(E,rho,I,A,kappa,...
                xlocs(elenodes(i,:)),options);




            % Master Stiffness
            for j = 1:n_kap
                K(dofs_in_ele,dofs_in_ele,j) = K(dofs_in_ele,dofs_in_ele,j) + ...
                    K_ele(:,:,j);
            end

            % Master Mass
            if strcmpi(options.elementType,'unconjugatedBlochOperator')
                for j = 1:n_kap
                    M(dofs_in_ele,dofs_in_ele,j) = M(dofs_in_ele,dofs_in_ele,j) + ...
                        M_ele(:,:,j);
                end
            else
                % Master Mass
                M(dofs_in_ele,dofs_in_ele) = M(dofs_in_ele,dofs_in_ele) + M_ele;
            end
        end

%         if strcmpi(options.elementType,'unconjugatedBlochOperator')
%             figure(99);clf
%             subplot(1,2,1);
%             bar3((K_ele(:,:,end)))
%             subplot(1,2,2);
%             bar3(imag(K_ele(:,:,end)))
%             pause
%         end




    % symbolic version of code
    case 'symbolic'

    % Default options
    defaults.elementType = 'standard'; 

    % fill in unspecified options fields with the default values
    options = setstructfields(defaults,options);

    % Mesh info
    Es              = Model_Info.Es;
    rhos            = Model_Info.rhos;
    Is              = Model_Info.Is;
    As              = Model_Info.As;
    pattern         = Model_Info.pattern;
    n_dof           = Model_Info.n_dof;
    elenodes        = Model_Info.elenodes;
    xlocs           = Model_Info.xlocs;

    n_kap  = size(kappa,2);
    n_eles = size(elenodes,1);

    K = zeros(n_dof,n_dof,n_kap);
    if strcmpi(options.elementType,'unconjugatedBlochOperator')
        M = zeros(n_dof,n_dof,n_kap);
    else
        M = zeros(n_dof,n_dof);
    end

    % symbolic element derivation
    syms x l k

    % 2 node Euler-Bernoulli Beam Shape Functions
    N = [1 - 3*x^2/l^2 + 2*x^3/l^3,...
         x - 2*x^2/l + x^3/l^2,...
         3*x^2/l^2 - 2*x^3/l^3,...
         -x^2/l + x^3/l^2]...
         *exp(1i*k*x);

    % 1st Derivative
    N_x = diff(N,x);

    % 2nd Derivative (B-matrix)
    B  = diff(N_x,x);

    % bloch-operator stiffness matrix
    switch options.elementType
        case 'BlochOperator'
            K_ele_sym = simplify(int(B'*B,x,0,l));
            M_ele_sym = simplify(int(N'*N,x,0,l));
        case 'unconjugatedBlochOperator'
            K_ele_sym = simplify(int(B.'*B,x,0,l));
            M_ele_sym = simplify(int(N.'*N,x,0,l));
        case 'standard'
            K_ele_sym = simplify(int(B'*B,x,0,l));
            M_ele_sym = simplify(int(N'*N,x,0,l));
    end

    for i = 1:n_eles
        i
        dofs_in_ele = zeros(4,1);
        for j = 1:length(elenodes(i,:))
            en = elenodes(i,j);
            dofs_in_ele(2*(j-1)+1:2*j) = [2*(en-1)+1:2*en];
        end

        rho = rhos{pattern(i)};
        E   = Es{pattern(i)};
        I   = Is{pattern(i)};
        A   = As{pattern(i)};

        % sub length into symbolic element mass and stiffness matrices
        ln = abs(xlocs(elenodes(i,2))-xlocs(elenodes(i,1)));
        K_ele_sym = eval(subs(K_ele_sym,l,ln));
        M_ele_sym = eval(subs(M_ele_sym,l,ln));

        % Master Stiffness
        for j = 1:n_kap
            K(dofs_in_ele,dofs_in_ele,j) = K(dofs_in_ele,dofs_in_ele,j) + ...
                E*I*eval(subs(K_ele_sym,k,kappa(j)));
        end

        % Master Mass
        if strcmpi(options.elementType,'unconjugatedBlochOperator')
            for j = 1:n_kap
                M(dofs_in_ele,dofs_in_ele,j) = M(dofs_in_ele,dofs_in_ele,j) + ...
                    rho*A*eval(subs(M_ele_sym,k,kappa(j)));
            end
        else
            % Master Mass
            M(dofs_in_ele,dofs_in_ele) = M(dofs_in_ele,dofs_in_ele) + rho*A*M_ele_sym;
        end
    end
end