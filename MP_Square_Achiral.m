clc
clear
format long

%% Input variables
theta=30;                                                                  % Waviness Angle
theta=deg2rad(theta);

%% Material properties
E=71e9;                                                                    % Young's modulus (Pa)
G=27e9;                                                                    % Shear modulus
v=0.32;                                                                    % Poisson's ratio
ro=2700;                                                                   % Density (kg/m^3)
L=0.015;                                                                   % Beam length (m)
t=0.001;                                                                   % Beam thickness
h=1;                                                                       % Beam height
Iz=h*t^3/12;                                                               % Moment of inertia (m^4)
A=t*h;                                                                     % Cross section area (m^2)

%% Element stiffness matrix (Euler-Bernoulli beam)
n=2^5;                                                                     % Number of elements for each ligament
if theta==0
    L_elm=L/n;                                                             % Element length
else
    L_elm=L*theta/n/sin(theta);
end
a=L_elm/2;
K_elm=[A*E/2/a,           0,             0, -A*E/2/a,             0,             0;
          0,  3*E*Iz/2/a^3,  3*E*Iz/2/a^2,        0, -3*E*Iz/2/a^3,  3*E*Iz/2/a^2;
          0,  3*E*Iz/2/a^2,      2*E*Iz/a,        0, -3*E*Iz/2/a^2,        E*Iz/a;
   -A*E/2/a,             0,             0,  A*E/2/a,             0,             0;
          0, -3*E*Iz/2/a^3, -3*E*Iz/2/a^2,        0,  3*E*Iz/2/a^3, -3*E*Iz/2/a^2;
          0,  3*E*Iz/2/a^2,        E*Iz/a,        0, -3*E*Iz/2/a^2,     2*E*Iz/a];

%% Stiffness Matrix of one ligament
% Element orientation
ori_elm=zeros(n,1);
for j=1:n
    ori_elm(j)=theta*(1-1/n-(j-1)*2/n);
end

% Stiffness Matrix of every Element
K_OA1=zeros(n*6);
for i=1:n
    T=[cos(ori_elm(i)), sin(ori_elm(i)), 0, 0, 0, 0;
        -sin(ori_elm(i)), cos(ori_elm(i)), 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, cos(ori_elm(i)), sin(ori_elm(i)), 0;
        0, 0, 0, -sin(ori_elm(i)), cos(ori_elm(i)), 0;
        0, 0, 0, 0, 0, 1];

    Ke=T'*K_elm*T;
    K_OA1(6*(i-1)+1:6*(i-1)+6,6*(i-1)+1:6*(i-1)+6)=Ke;
end

% Assembling matrix
I=[1,0,0;0,1,0;0,0,1];                                                     
As=zeros(n*2*3,(n+1)*3);
As(1:3,1:3)=I;
for i=2:n*2
    j=floor(i/2);
    As((i-1)*3+1:(i-1)*3+3,j*3+1:j*3+3)=I;
end
K_OA2=As'*K_OA1*As;

% Simplified Stiffness Matrix of one ligament
K1(:,1:3)=K_OA2(:,1:3); K1(:,4:6)=K_OA2(:,3*n+1:3*n+3); K1(:,7:3*n+3)=K_OA2(:,4:3*n);
K2(1:3,:)=K1(1:3,:);K2(4:6,:)=K1(3*n+1:3*n+3,:);K2(7:3*n+3,:)=K1(4:3*n,:);
aa=K2(1:6,1:6);
ab=K2(1:6,7:end);
ba=K2(7:end,1:6);
bb=K2(7:end,7:end);
K_OA=aa-ab/bb*ba;
K_OA(abs(K_OA)<1e-3)=0;

%% Global stiffness matrix K_OAtoL
% Connectivity matrix
cnct=[2,10; 6,10; 1,10; 9,10; 6,11; 3,11; 7,11; 1,11; 1,12; 7,12; 4,12; 8,12; 9,13; 1,13; 8,13; 5,13;];
% Beam orientation
ori_lgm=[pi, 3*pi/2, 0, pi/2, pi, 3*pi/2, 0, pi/2, pi, 3*pi/2, 0, pi/2, pi, 3*pi/2, 0, pi/2];
K_OAtoL=zeros(13*3);
for i=1:length(ori_lgm)
    T=[cos(ori_lgm(i)), sin(ori_lgm(i)), 0, 0, 0, 0;
        -sin(ori_lgm(i)), cos(ori_lgm(i)), 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, cos(ori_lgm(i)), sin(ori_lgm(i)), 0;
        0, 0, 0, -sin(ori_lgm(i)), cos(ori_lgm(i)), 0;
        0, 0, 0, 0, 0, 1];

    Ke=T'*K_OA*T;
    Ke(abs(Ke)<1e-3)=0;
    node1=cnct(i,1);
    node2=cnct(i,2);
    K_OAtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node1-1)+1:3*(node1-1)+3)=...
    K_OAtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node1-1)+1:3*(node1-1)+3)+Ke(1:3,1:3);
    K_OAtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node2-1)+1:3*(node2-1)+3)=...
    K_OAtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node2-1)+1:3*(node2-1)+3)+Ke(1:3,4:6);
    K_OAtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node1-1)+1:3*(node1-1)+3)=...
    K_OAtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node1-1)+1:3*(node1-1)+3)+Ke(4:6,1:3);
    K_OAtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node2-1)+1:3*(node2-1)+3)=...
    K_OAtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node2-1)+1:3*(node2-1)+3)+Ke(4:6,4:6);

end

%% Global stiffness matrix K_OAtoH
nS=8;                                                                      % Number of joints on the surface of RVE
aa=K_OAtoL(1:3*nS+3,1:3*nS+3);
ab=K_OAtoL(1:3*nS+3,3*nS+4:end);
ba=K_OAtoL(3*nS+4:end,1:3*nS+3);
bb=K_OAtoL(3*nS+4:end,3*nS+4:end);
K_OAtoH=aa-ab/bb*ba;
K_OAtoH(abs(K_OAtoH)<1e-3)=0;

%% Taylor's Expansion
L1=2*L;
L2=sqrt(2)*L;
ori1=[0,pi/2,pi,3*pi/2];
ori2=[pi/4,3*pi/4,5*pi/4,7*pi/4];
syms e11 e22 e12 e21 k13 k23 phi0 u0x u0y real
A12=(e12-e21)/2;
u12=e21-phi0;
u21=e12+phi0;
psi0=A12;
psi1=A12-(k13*L1*cos(ori1(1))+k23*L1*sin(ori1(1)));
psi2=A12-(k13*L1*cos(ori1(2))+k23*L1*sin(ori1(2)));
psi3=A12-(k13*L1*cos(ori1(3))+k23*L1*sin(ori1(3)));
psi4=A12-(k13*L1*cos(ori1(4))+k23*L1*sin(ori1(4)));
u1x=u0x+e11*L1*cos(ori1(1))+u12*L1*sin(ori1(1));
u1y=u0y+e22*L1*sin(ori1(1))+u21*L1*cos(ori1(1));
u2x=u0x+e11*L1*cos(ori1(2))+u12*L1*sin(ori1(2));
u2y=u0y+e22*L1*sin(ori1(2))+u21*L1*cos(ori1(2));
u3x=u0x+e11*L1*cos(ori1(3))+u12*L1*sin(ori1(3));
u3y=u0y+e22*L1*sin(ori1(3))+u21*L1*cos(ori1(3));
u4x=u0x+e11*L1*cos(ori1(4))+u12*L1*sin(ori1(4));
u4y=u0y+e22*L1*sin(ori1(4))+u21*L1*cos(ori1(4));

psi5=A12-(k13*L2*cos(ori2(1))+k23*L2*sin(ori2(1)));
psi6=A12-(k13*L2*cos(ori2(2))+k23*L2*sin(ori2(2)));
psi7=A12-(k13*L2*cos(ori2(3))+k23*L2*sin(ori2(3)));
psi8=A12-(k13*L2*cos(ori2(4))+k23*L2*sin(ori2(4)));
u5x=u0x+e11*L2*cos(ori2(1))+u12*L2*sin(ori2(1));
u5y=u0y+e22*L2*sin(ori2(1))+u21*L2*cos(ori2(1));
u6x=u0x+e11*L2*cos(ori2(2))+u12*L2*sin(ori2(2));
u6y=u0y+e22*L2*sin(ori2(2))+u21*L2*cos(ori2(2));
u7x=u0x+e11*L2*cos(ori2(3))+u12*L2*sin(ori2(3));
u7y=u0y+e22*L2*sin(ori2(3))+u21*L2*cos(ori2(3));
u8x=u0x+e11*L2*cos(ori2(4))+u12*L2*sin(ori2(4));
u8y=u0y+e22*L2*sin(ori2(4))+u21*L2*cos(ori2(4));

% Solve for u0x, u0y, psi0
Kuu=K_OAtoH(1:3,1:3); Kuk=K_OAtoH(1:3,4:end);
Kp=-Kuu\Kuk;
Kp(abs(Kp)<1e-3)=0;

dk=[u1x;u1y;psi1;u2x;u2y;psi2;u3x;u3y;psi3;u4x;u4y;psi4;u5x;u5y;psi5;u6x;u6y;psi6;u7x;u7y;psi7;u8x;u8y;psi8;];
du=[u0x;u0y;psi0];
s=solve(Kp*dk==du,u0x,u0y,phi0);
phi0s=vpa(s.phi0,4);
u0xs=vpa(s.u0x,4);
u0ys=vpa(s.u0y,4);

p0=subs(psi0,phi0,phi0s);
p1=subs(psi1,phi0,phi0s);
p2=subs(psi2,phi0,phi0s);
p3=subs(psi3,phi0,phi0s);
p4=subs(psi4,phi0,phi0s);
p5=subs(psi5,phi0,phi0s);
p6=subs(psi6,phi0,phi0s);
p7=subs(psi7,phi0,phi0s);
p8=subs(psi8,phi0,phi0s);
u1x1=subs(u1x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u1y1=subs(u1y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u2x1=subs(u2x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u2y1=subs(u2y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u3x1=subs(u3x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u3y1=subs(u3y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u4x1=subs(u4x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u4y1=subs(u4y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u5x1=subs(u5x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u5y1=subs(u5y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u6x1=subs(u6x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u6y1=subs(u6y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u7x1=subs(u7x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u7y1=subs(u7y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u8x1=subs(u8x,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);
u8y1=subs(u8y,[u0x,u0y,phi0],[s.u0x,s.u0y,phi0s]);

%% Micropolar Elasticity Tensor
d=[u0xs;u0ys;p0;u1x1;u1y1;p1;u2x1;u2y1;p2;u3x1;u3y1;p3;u4x1;u4y1;p4;u5x1;u5y1;p5;u6x1;u6y1;p6;u7x1;u7y1;p7;u8x1;u8y1;p8;];
V=8*L^2;                                                                   % Volume of RVE
w=d'*K_OAtoH*d/(2*V);                                                      % Strain Energy Density
c1=diff(w,e11,e11);c2=diff(w,e11,e22);c3=diff(w,e11,e12);
c4=diff(w,e11,e21);c5=diff(w,e11,k13);c6=diff(w,e11,k23);
c7=diff(w,e22,e22);c8=diff(w,e22,e12);c9=diff(w,e22,e21);
c10=diff(w,e22,k13);c11=diff(w,e22,k23);c12=diff(w,e12,e12);
c13=diff(w,e12,e21);c14=diff(w,e12,k13);c15=diff(w,e12,k23);
c16=diff(w,e21,e21);c17=diff(w,e21,k13);c18=diff(w,e21,k23);
c19=diff(w,k13,k13);c20=diff(w,k13,k23);c21=diff(w,k23,k23);
Q=[c1, c2, c3, c4, c5, c6;
   c2, c7, c8, c9,c10,c11;
   c3, c8,c12,c13,c14,c15;
   c4, c9,c13,c16,c17,c18;
   c5,c10,c14,c17,c19,c20;
   c6,c11,c15,c18,c20,c21];
Q=double(Q);
Q(abs(Q)<1e-3)=0

%% Mechanical Properties
C=Q(1:4,1:4);
D=Q(5:6,5:6);
B=Q(1:4,5:6);
Qi=inv(Q);
Em=1/Qi(1,1);                                                              % Effective Young's Modulus
vm=-Qi(2,1)/Qi(1,1);                                                       % Effective Poisson's Ratio
Gm=1/(Qi(3,3)+Qi(3,4))/2;                                                  % Effective Shear Modulus
Bm=1/(Qi(1,1)+Qi(1,2))/2;                                                  % Effective Bulk Modulus
yita=(Qi(1,3)+Qi(1,4))/Qi(1,1);                                            % Axial-Shear Coupling Coefficient
