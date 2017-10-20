function Me = Tet10Mass(xyz,rho)

% GQ integration rule
rule=14;

% preallocate stiffness array
Me = zeros(30,30);
 
 for k = 1:abs(rule)
     
     % get GQ coordinate and weight
     [zeta,wk] = TetGaussRuleInfo(rule,k);
     
     % get shape function derivatives and Jacobian determinant
     [N,Jdet]=Tet10ShapeFun(xyz,zeta);
     
     if Jdet<=0, 
         disp('Tet10Stiffness: Neg ');
     end
     
     % form element strain-displacement matrix, B
     thr = (0:3:27);
     Nmat = zeros(3,30);
     Nmat(1,(thr+1)) = N;
     Nmat(2,(thr+2)) = N;
     Nmat(3,(thr+3)) = N;
     
     Me = Me + (wk*Jdet/6)*rho*Nmat'*Nmat;
 end