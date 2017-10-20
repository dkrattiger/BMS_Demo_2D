function Me = Tet4Mass(xyz,rho)

% xyz = [0,0,0;...
%        1,0,0;...
%        0,1,0;...
%        0,0,1];

% degrees of freedom per node
n_dpn = 3;

% extract coordinates
x1 = xyz(1,1);
x2 = xyz(2,1);
x3 = xyz(3,1);
x4 = xyz(4,1);

y1 = xyz(1,2);
y2 = xyz(2,2);
y3 = xyz(3,2);
y4 = xyz(4,2);

z1 = xyz(1,3);
z2 = xyz(2,3);
z3 = xyz(3,3);
z4 = xyz(4,3);

J = [ones(1,4);xyz'];

V01 = (1/6)*(x2*(y3*z4-y4*z3) + x3*(y4*z2-y2*z4) + x4*(y2*z3-y3*z2));
V02 = (1/6)*(x1*(y4*z3-y3*z4) + x3*(y1*z4-y4*z1) + x4*(y3*z1-y1*z3));
V03 = (1/6)*(x1*(y2*z4-y4*z2) + x2*(y4*z1-y1*z4) + x4*(y1*z2-y2*z1));
V04 = (1/6)*(x1*(y3*z2-y2*z3) + x2*(y1*z3-y3*z1) + x3*(y2*z1-y1*z2));
V = V01+V02+V03+V04

a1 = (y4-y2)*(z3-z2) - (y3-y2)*(z4-z2);
a2 = (y3-y1)*(z4-z3) - (y3-y4)*(z1-z3);
a3 = (y2-y4)*(z1-z4) - (y1-y4)*(z2-z4);
a4 = (y1-y3)*(z2-z1) - (y1-y2)*(z3-z1);

b1 = (x3-x2)*(z4-z2) - (x4-x2)*(z3-z2);
b2 = (x4-x3)*(z3-z1) - (x1-x3)*(z3-z4);
b3 = (x1-x4)*(z2-z4) - (x2-x4)*(z1-z4);
b4 = (x2-x1)*(z1-z3) - (x3-x1)*(z1-z2);

c1 = (x4-x2)*(y3-y2) - (x3-x2)*(y4-y2);
c2 = (x3-x1)*(y4-y3) - (x3-x4)*(y1-y3);
c3 = (x2-x4)*(y1-y4) - (x1-x4)*(y2-y4);
c4 = (x1-x3)*(y2-y1) - (x1-x2)*(y3-y1);

Rho1 = inv(J)
