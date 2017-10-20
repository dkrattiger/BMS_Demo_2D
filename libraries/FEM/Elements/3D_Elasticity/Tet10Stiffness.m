function Ke = Tet10Stiffness(xyz,Emat)

% GQ integration rule
rule=4;

% preallocate stiffness array
Ke = zeros(30,30);
 
 for k = 1:abs(rule)
     
     % get GQ coordinate and weight
     [zeta,wk] = TetGaussRuleInfo(rule,k);
     
     % get shape function derivatives and Jacobian determinant
     [Nx,Ny,Nz,Jdet]=Tet10ShapeFunDer(xyz,zeta);
     
     if Jdet<=0, 
         disp('Tet10Stiffness: Neg ');
     end
     
     % form element strain-displacement matrix, B
     thr = (0:3:27);
     B = zeros(6,30);
     B(1,(thr+1)) = Nx;
     B(2,(thr+2)) = Ny;
     B(3,(thr+3)) = Nz;
     B(4,(thr+1)) = Ny;
     B(4,(thr+2)) = Nx;
     B(5,(thr+2)) = Nz;
     B(5,(thr+3)) = Ny;
     B(6,(thr+1)) = Nz;
     B(6,(thr+3)) = Nx;
     
     Ke = Ke + (wk/(6*Jdet))*B'*Emat*B;
 end