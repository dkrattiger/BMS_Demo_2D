function [N,Jdet] = Tet10ShapeFun(xyz,zeta)

x1 = xyz(1,1);
x2 = xyz(2,1);
x3 = xyz(3,1);
x4 = xyz(4,1);
x5 = xyz(5,1);
x6 = xyz(6,1);
x7 = xyz(7,1);
x8 = xyz(8,1);
x9 = xyz(9,1);
x10 = xyz(10,1);

y1 = xyz(1,2);
y2 = xyz(2,2);
y3 = xyz(3,2);
y4 = xyz(4,2);
y5 = xyz(5,2);
y6 = xyz(6,2);
y7 = xyz(7,2);
y8 = xyz(8,2);
y9 = xyz(9,2);
y10 = xyz(10,2);

z1 = xyz(1,3);
z2 = xyz(2,3);
z3 = xyz(3,3);
z4 = xyz(4,3);
z5 = xyz(5,3);
z6 = xyz(6,3);
z7 = xyz(7,3);
z8 = xyz(8,3);
z9 = xyz(9,3);
z10 = xyz(10,3);

zeta1 = zeta(1);
zeta2 = zeta(2);
zeta3 = zeta(3);
zeta4 = zeta(4);


% Define Jacobian terms
Jx1=4*(x1*(zeta1-1/4)+x5*zeta2+x7*zeta3+x8*zeta4);
Jy1=4*(y1*(zeta1-1/4)+y5*zeta2+y7*zeta3+y8*zeta4);
Jz1=4*(z1*(zeta1-1/4)+z5*zeta2+z7*zeta3+z8*zeta4);
Jx2=4*(x5*zeta1+x2*(zeta2-1/4)+x6*zeta3+x9*zeta4);
Jy2=4*(y5*zeta1+y2*(zeta2-1/4)+y6*zeta3+y9*zeta4);
Jz2=4*(z5*zeta1+z2*(zeta2-1/4)+z6*zeta3+z9*zeta4);
Jx3=4*(x7*zeta1+x6*zeta2+x3*(zeta3-1/4)+x10*zeta4);
Jy3=4*(y7*zeta1+y6*zeta2+y3*(zeta3-1/4)+y10*zeta4);
Jz3=4*(z7*zeta1+z6*zeta2+z3*(zeta3-1/4)+z10*zeta4);
Jx4=4*(x8*zeta1+x9*zeta2+x10*zeta3+x4*(zeta4-1/4));
Jy4=4*(y8*zeta1+y9*zeta2+y10*zeta3+y4*(zeta4-1/4));
Jz4=4*(z8*zeta1+z9*zeta2+z10*zeta3+z4*(zeta4-1/4));
Jx12=Jx1-Jx2; Jx13=Jx1-Jx3; Jx14=Jx1-Jx4; Jx23=Jx2-Jx3;
Jx24=Jx2-Jx4; Jx34=Jx3-Jx4; Jy12=Jy1-Jy2; Jy13=Jy1-Jy3;
Jy14=Jy1-Jy4; Jy23=Jy2-Jy3; Jy24=Jy2-Jy4; Jy34=Jy3-Jy4;
Jz12=Jz1-Jz2; Jz13=Jz1-Jz3; Jz14=Jz1-Jz4;
Jz23=Jz2-Jz3; Jz24=Jz2-Jz4; Jz34=Jz3-Jz4;
Jx21=-Jx12; Jx31=-Jx13; Jx41=-Jx14; Jx32=-Jx23; Jx42=-Jx24;
Jx43=-Jx34; Jy21=-Jy12; Jy31=-Jy13; Jy41=-Jy14; Jy32=-Jy23;
Jy42=-Jy24; Jy43=-Jy34; Jz21=-Jz12; Jz31=-Jz13; Jz41=-Jz14;
Jz32=-Jz23; Jz42=-Jz24; Jz43=-Jz34;

% Jacobian determinant
Jdet=Jx21*(Jy23*Jz34-Jy34*Jz23)+Jx32*(Jy34*Jz12-Jy12*Jz34)+...
           Jx43*(Jy12*Jz23-Jy23*Jz12);

% Jacobian inverse terms
a1=Jy42*Jz32-Jy32*Jz42; a2=Jy31*Jz43-Jy34*Jz13;
a3=Jy24*Jz14-Jy14*Jz24; a4=Jy13*Jz21-Jy12*Jz31;
b1=Jx32*Jz42-Jx42*Jz32; b2=Jx43*Jz31-Jx13*Jz34;
b3=Jx14*Jz24-Jx24*Jz14; b4=Jx21*Jz13-Jx31*Jz12;
c1=Jx42*Jy32-Jx32*Jy42; c2=Jx31*Jy43-Jx34*Jy13;
c3=Jx24*Jy14-Jx14*Jy24; c4=Jx13*Jy21-Jx12*Jy31;
 
N = [zeta1*(2*zeta1-1),...
     zeta2*(2*zeta2-1),...
     zeta3*(2*zeta3-1),...
     zeta4*(2*zeta4-1),...
     4*zeta1*zeta2,...
     4*zeta2*zeta3,...
     4*zeta3*zeta1,...
     4*zeta1*zeta4,...
     4*zeta2*zeta4,...
     4*zeta3*zeta4];