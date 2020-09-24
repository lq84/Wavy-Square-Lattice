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
for i=1:n/2
    ori_elm(i)=theta*(1-2/n-(i-1)*4/n);
    ori_elm(i+n/2)=-theta*(1-2/n-(i-1)*4/n);
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
    K_OA1(6*(i-1)+1:6*(i-1)+6,6*(i-1)+1:6*(i-1)+6)=T'*K_elm*T;
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
K1(:,1:3)=K_OA2(:,1:3);K1(:,4:6)=K_OA2(:,3*n+1:3*n+3);K1(:,7:3*n+3)=K_OA2(:,4:3*n);
K2(1:3,:)=K1(1:3,:);K2(4:6,:)=K1(3*n+1:3*n+3,:);K2(7:3*n+3,:)=K1(4:3*n,:);
aa=K2(1:6,1:6);
ab=K2(1:6,7:3*n+3);
ba=K2(7:3*n+3,1:6);
bb=K2(7:3*n+3,7:3*n+3);
K_OA=aa-ab*inv(bb)*ba;
K_OA(abs(K_OA)<1e-3)=0;                                                      

%% Global stiffness matrix
% Stiffness Matrix of four ligaments
ori_lgm=[0,pi/2,pi,3/2*pi];
K_OABCD1=zeros(4*6);
for i=1:4
    T=[cos(ori_lgm(i)), sin(ori_lgm(i)), 0, 0, 0, 0;
    -sin(ori_lgm(i)), cos(ori_lgm(i)), 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, cos(ori_lgm(i)), sin(ori_lgm(i)), 0;
    0, 0, 0, -sin(ori_lgm(i)), cos(ori_lgm(i)), 0;
    0, 0, 0, 0, 0, 1];
    Ke=T'*K_OA*T;
    Ke(abs(Ke)<1e-3)=0;
    K_OABCD1(6*(i-1)+1:6*(i-1)+6,6*(i-1)+1:6*(i-1)+6)=Ke;
end

% Assembling matrix
I=[1,0,0;0,1,0;0,0,1];
As=zeros(4*2*3,(4+1)*3);
for i=1:4
    As(((i-1)*2)*3+1:((i-1)*2)*3+3,1:3)=I;
    j=2;
    jj=floor(j/2);
    As((j-1+(i-1)*2)*3+1:(j-1+(i-1)*2)*3+3,...
        (jj+(i-1))*3+1:(jj+(i-1))*3+3)=I;
end
K_OABCD=As'*K_OABCD1*As;

%% Taylor's Expansion
L_lgm=L;
ori_lgm=[0,pi/2,pi,3/2*pi];
syms e11 e22 e12 e21 k13 k23 phi0 u0x u0y real
A12=(e12-e21)/2;
u12=e21-phi0;
u21=e12+phi0;
psi0=A12;
psi1=A12-k13*L_lgm;
psi2=A12-k23*L_lgm;
psi3=A12+k13*L_lgm;
psi4=A12+k23*L_lgm;

u1x=u0x+e11*L_lgm*cos(ori_lgm(1))+u12*L_lgm*sin(ori_lgm(1));
u1y=u0y+e22*L_lgm*sin(ori_lgm(1))+u21*L_lgm*cos(ori_lgm(1));
u2x=u0x+e11*L_lgm*cos(ori_lgm(2))+u12*L_lgm*sin(ori_lgm(2));
u2y=u0y+e22*L_lgm*sin(ori_lgm(2))+u21*L_lgm*cos(ori_lgm(2));
u3x=u0x+e11*L_lgm*cos(ori_lgm(3))+u12*L_lgm*sin(ori_lgm(3));
u3y=u0y+e22*L_lgm*sin(ori_lgm(3))+u21*L_lgm*cos(ori_lgm(3));
u4x=u0x+e11*L_lgm*cos(ori_lgm(4))+u12*L_lgm*sin(ori_lgm(4));
u4y=u0y+e22*L_lgm*sin(ori_lgm(4))+u21*L_lgm*cos(ori_lgm(4));

% Solve for u0x, u0y, psi0
Kuu=K_OABCD(1:3,1:3);Kuk=K_OABCD(1:3,4:15);
Kp=-Kuu\Kuk;
Kp(abs(Kp)<1e-3)=0;

du=[u0x;u0y;psi0];
dk=[u1x;u1y;psi1;u2x;u2y;psi2;u3x;u3y;psi3;u4x;u4y;psi4];
s=solve(Kp*dk==du,u0x,u0y,phi0);
sw0=vpa(s.phi0,4);
u0x1=s.u0x;
u0y1=s.u0y;

u1x1=subs(u1x,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u1y1=subs(u1y,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u2x1=subs(u2x,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u2y1=subs(u2y,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u3x1=subs(u3x,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u3y1=subs(u3y,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u4x1=subs(u4x,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
u4y1=subs(u4y,[u0x,u0y,phi0],[u0x1,u0y1,sw0]);
p0=subs(psi0,phi0,sw0);
p1=subs(psi1,phi0,sw0);
p2=subs(psi2,phi0,sw0);
p3=subs(psi3,phi0,sw0);
p4=subs(psi4,phi0,sw0);

%% Micropolar Elasticity Tensor
d=[u0x1;u0y1;p0;u1x1;u1y1;p1;u2x1;u2y1;p2;u3x1;u3y1;p3;u4x1;u4y1;p4];
w=d'*K_OABCD*d/L^2/4;                                                      % Strain Energy Density
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

