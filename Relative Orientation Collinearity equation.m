
% Relative Orientation based on the collinearity equations of terrestrial Photogrammetry.
% input:
% - image coordinates [mm] :xx,yy [in p.p. system](here 6 points, common in both images)
% - c: focal length [mm]                                                
%-initial exterior orienattion parameters                       
% output:                                                               
% - adjusted relative orienattion parameters (by, bz, omega, phi, kappa)                           


clear all
format long 
%left photo coordinates(mm)
left = [9.44 33.25
       23.85 32.21
       36.46 33.20
       12.66 -24.04
       20.70 -7.58
       43.99 -23.80];
%right photo coordinates(mm)
right = [-29.40 31.44
        -15.24 30.34
        -2.64 31.32
        -32.93 -26.05
        -22.38 -9.54
        -1.82 -25.84];
%focal length
c = 99.91;
%the approximate exterior orientation parameters of the left image
X0l = 0;
Y0l = 0;
Z0l = 0;
omega_l = 0;
phi_l = 0;
kappa_l = 0;
%the approximate exterior orientation parameters of the right image
bx = c;
by = 0;
bz = 0;
omega = 0;
phi = 0;
kappa = 0;
%approximate reference point coordinates
X = zeros(6,1);
for i = 1:6
    %parallax
    px = left(i,1) - right(i,1);
    X(i,1) = left(i,1)*(bx/px);
end
Y = zeros(6,1);
for i = 1:6
    px = left(i,1) - right(i,1);
    Y(i,1) = left(i,2)*(bx/px);
end
Z = zeros(6,1);
H = zeros(6,1);
for i =1:6
    px = left(i,1) - right(i,1);
    H(i) = (c*bx)/px;
    Z(i,1) = Z0l - H(i);
end
geo = zeros(6,3);
geo(1:6,1) = X;
geo(1:6,2) = Y;
geo(1:6,3) = Z;


omega = omega*pi()/200;
phi = phi*pi()/200;
kappa = kappa*pi()/200;
% ñ^cc
r = 636620;
% Rùöê
R = zeros(3,3);
R(1,1) = cos(phi)*cos(kappa);
R(1,2) = sin(omega)*sin(phi)*cos(kappa)+cos(omega)*sin(kappa);
R(1,3) = -cos(omega)*sin(phi)*cos(kappa)+sin(omega)*sin(kappa);
R(2,1) = -cos(phi)*sin(kappa);
R(2,2) = -sin(omega)*sin(phi)*sin(kappa)+cos(omega)*cos(kappa);
R(2,3) = cos(omega)*sin(phi)*sin(kappa)+sin(omega)*cos(kappa);
R(3,1) = sin(phi);
R(3,2) = -sin(omega)*cos(phi);
R(3,3) = cos(omega)*cos(phi);
%partial derivative Rùöê(ö)
f = zeros(3,3);
f(1,1) = -sin(phi)*cos(kappa);
f(1,2) = sin(omega)*cos(phi)*cos(kappa);
f(1,3) = -cos(omega)*cos(phi)*cos(kappa);
f(2,1) = sin(phi)*sin(kappa);
f(2,2) = -sin(omega)*cos(phi)*sin(kappa);
f(2,3) = cos(omega)*cos(phi)*sin(kappa);
f(3,1) = cos(phi);
f(3,2) = sin(omega)*sin(phi);
f(3,3) = -cos(omega)*sin(phi);
%numerator for x in the collinearity equation for every ground control point (left image)
A1_l = zeros(6,1);
A1_l(1,1) = R(1,1)*(geo(1,1)-X0l)+R(1,2)*(geo(1,2)-Y0l)+R(1,3)*(geo(1,3)-Z0l);
A1_l(2,1) = R(1,1)*(geo(2,1)-X0l)+R(1,2)*(geo(2,2)-Y0l)+R(1,3)*(geo(2,3)-Z0l);
A1_l(3,1) = R(1,1)*(geo(3,1)-X0l)+R(1,2)*(geo(3,2)-Y0l)+R(1,3)*(geo(3,3)-Z0l);
A1_l(4,1) = R(1,1)*(geo(4,1)-X0l)+R(1,2)*(geo(4,2)-Y0l)+R(1,3)*(geo(4,3)-Z0l);
A1_l(5,1) = R(1,1)*(geo(5,1)-X0l)+R(1,2)*(geo(5,2)-Y0l)+R(1,3)*(geo(5,3)-Z0l);
A1_l(6,1) = R(1,1)*(geo(6,1)-X0l)+R(1,2)*(geo(6,2)-Y0l)+R(1,3)*(geo(6,3)-Z0l);
%numerator for y in the collinearity equation for every ground control point (left image)
A2_l = zeros(6,1);
A2_l(1,1) = R(2,1)*(geo(1,1)-X0l)+R(2,2)*(geo(1,2)-Y0l)+R(2,3)*(geo(1,3)-Z0l);
A2_l(2,1) = R(2,1)*(geo(2,1)-X0l)+R(2,2)*(geo(2,2)-Y0l)+R(2,3)*(geo(2,3)-Z0l);
A2_l(3,1) = R(2,1)*(geo(3,1)-X0l)+R(2,2)*(geo(3,2)-Y0l)+R(2,3)*(geo(3,3)-Z0l);
A2_l(4,1) = R(2,1)*(geo(4,1)-X0l)+R(2,2)*(geo(4,2)-Y0l)+R(2,3)*(geo(4,3)-Z0l);
A2_l(5,1) = R(2,1)*(geo(5,1)-X0l)+R(2,2)*(geo(5,2)-Y0l)+R(2,3)*(geo(5,3)-Z0l);
A2_l(6,1) = R(2,1)*(geo(6,1)-X0l)+R(2,2)*(geo(6,2)-Y0l)+R(2,3)*(geo(6,3)-Z0l);
%denominator in the collinearity equation for every ground control point (left image)
P_l = zeros(6,1);
P_l(1,1) = R(3,1)*(geo(1,1)-X0l)+R(3,2)*(geo(1,2)-Y0l)+R(3,3)*(geo(1,3)-Z0l);
P_l(2,1) = R(3,1)*(geo(2,1)-X0l)+R(3,2)*(geo(2,2)-Y0l)+R(3,3)*(geo(2,3)-Z0l);
P_l(3,1) = R(3,1)*(geo(3,1)-X0l)+R(3,2)*(geo(3,2)-Y0l)+R(3,3)*(geo(3,3)-Z0l);
P_l(4,1) = R(3,1)*(geo(4,1)-X0l)+R(3,2)*(geo(4,2)-Y0l)+R(3,3)*(geo(4,3)-Z0l);
P_l(5,1) = R(3,1)*(geo(5,1)-X0l)+R(3,2)*(geo(5,2)-Y0l)+R(3,3)*(geo(5,3)-Z0l);
P_l(6,1) = R(3,1)*(geo(6,1)-X0l)+R(3,2)*(geo(6,2)-Y0l)+R(3,3)*(geo(6,3)-Z0l);
%numerator for x in the collinearity equation for every ground control point (right image)
A1_r = zeros(6,1);
A1_r(1,1) = R(1,1)*(geo(1,1)-bx)+R(1,2)*(geo(1,2)-by)+R(1,3)*(geo(1,3)-bz);
A1_r(2,1) = R(1,1)*(geo(2,1)-bx)+R(1,2)*(geo(2,2)-by)+R(1,3)*(geo(2,3)-bz);
A1_r(3,1) = R(1,1)*(geo(3,1)-bx)+R(1,2)*(geo(3,2)-by)+R(1,3)*(geo(3,3)-bz);
A1_r(4,1) = R(1,1)*(geo(4,1)-bx)+R(1,2)*(geo(4,2)-by)+R(1,3)*(geo(4,3)-bz);
A1_r(5,1) = R(1,1)*(geo(5,1)-bx)+R(1,2)*(geo(5,2)-by)+R(1,3)*(geo(5,3)-bz);
A1_r(6,1) = R(1,1)*(geo(6,1)-bx)+R(1,2)*(geo(6,2)-by)+R(1,3)*(geo(6,3)-bz);
%numerator for y in the collinearity equation for every ground control point (right image)

A2_r = zeros(6,1);
A2_r(1,1) = R(2,1)*(geo(1,1)-bx)+R(2,2)*(geo(1,2)-by)+R(2,3)*(geo(1,3)-bz);
A2_r(2,1) = R(2,1)*(geo(2,1)-bx)+R(2,2)*(geo(2,2)-by)+R(2,3)*(geo(2,3)-bz);
A2_r(3,1) = R(2,1)*(geo(3,1)-bx)+R(2,2)*(geo(3,2)-by)+R(2,3)*(geo(3,3)-bz);
A2_r(4,1) = R(2,1)*(geo(4,1)-bx)+R(2,2)*(geo(4,2)-by)+R(2,3)*(geo(4,3)-bz);
A2_r(5,1) = R(2,1)*(geo(5,1)-bx)+R(2,2)*(geo(5,2)-by)+R(2,3)*(geo(5,3)-bz);
A2_r(6,1) = R(2,1)*(geo(6,1)-bx)+R(2,2)*(geo(6,2)-by)+R(2,3)*(geo(6,3)-bz);
%denominator in the collinearity equation for every ground control point (right image)
P_r = zeros(6,1);
P_r(1,1) = R(3,1)*(geo(1,1)-bx)+R(3,2)*(geo(1,2)-by)+R(3,3)*(geo(1,3)-bz);
P_r(2,1) = R(3,1)*(geo(2,1)-bx)+R(3,2)*(geo(2,2)-by)+R(3,3)*(geo(2,3)-bz);
P_r(3,1) = R(3,1)*(geo(3,1)-bx)+R(3,2)*(geo(3,2)-by)+R(3,3)*(geo(3,3)-bz);
P_r(4,1) = R(3,1)*(geo(4,1)-bx)+R(3,2)*(geo(4,2)-by)+R(3,3)*(geo(4,3)-bz);
P_r(5,1) = R(3,1)*(geo(5,1)-bx)+R(3,2)*(geo(5,2)-by)+R(3,3)*(geo(5,3)-bz);
P_r(6,1) = R(3,1)*(geo(6,1)-bx)+R(3,2)*(geo(6,2)-by)+R(3,3)*(geo(6,3)-bz);
%matrix Á
A = zeros(24,23);

%partial derivatives with respect to x1 (left image)
A(1,6) = -c*(R(1,1)/P_l(1,1));
A(1,8) = c*(R(3,3)*A1_l(1,1)/P_l(1,1)^2);
%partial derivatives with respect to y1 (left image)
A(2,7) = -c*(R(2,2)/P_l(1,1));
A(2,8) = c*(R(3,3)*A2_l(1,1)/P_l(1,1)^2);
%partial derivatives with respect to x1 (right image)
A(3,2) = -c*(R(3,3)*A1_r(1,1)/P_r(1,1)^2);
A(3,3) = -c*(R(1,1)*A1_r(1,1)*A2_r(1,1)/P_r(1,1)^2);
A(3,4) = -c*((sin(kappa)*(A1_r(1,1)*A2_r(1,1)/P_r(1,1)^2))-cos(kappa)*(1+(A1_r(1,1)^2)/P_r(1,1)^2));
A(3,5) = -c*(A2_r(1,1)/P_r(1,1));
A(3,6) = -c*(R(1,1)/P_r(1,1));
A(3,8) = c*(R(3,3)*A1_r(1,1)/P_r(1,1)^2);
%partial derivatives with respect to y1 (right image)
A(4,1) = c*(R(2,2)/P_r(1,1));
A(4,2) = -c*(R(3,3)*A2_r(1,1)/P_r(1,1)^2);
A(4,3) = -(c/P_r(1,1)^2)*((R(1,1)*P_r(1,1)^2)+(R(1,1)*A2_r(1,1)^2));
A(4,4) = -c*(sin(kappa)*(1+(A2_r(1,1)^2)/P_r(1,1)^2)-cos(kappa)*A1_r(1,1)*A2_r(1,1)/P_r(1,1)^2);
A(4,5) = c*A1_r(1,1)/P_r(1,1);
A(4,7) = -c*R(2,2)/P_r(1,1);
A(4,8) = c*(R(3,3)*A2_r(1,1)/P_r(1,1)^2);
%partial derivatives with respect to x2 (left image)
A(5,9) = -c*(R(1,1)/P_l(2,1));
A(5,11) = c*(R(3,3)*A1_l(2,1)/P_l(2,1)^2);
%partial derivatives with respect to y2 (left image)
A(6,10) = -c*(R(2,2)/P_l(2,1));
A(6,11) = c*(R(3,3)*A2_l(2,1)/P_l(2,1)^2);
%partial derivatives with respect to x2 (right image)
A(7,2) = -c*(R(3,3)*A1_r(2,1)/P_r(2,1)^2);
A(7,3) = -c*(R(1,1)*A1_r(2,1)*A2_r(2,1)/P_r(2,1)^2);
A(7,4) = -c*((sin(kappa)*(A1_r(2,1)*A2_r(2,1)/P_r(2,1)^2))-cos(kappa)*(1+(A1_r(2,1)^2)/P_r(2,1)^2));
A(7,5) = -c*(A2_r(2,1)/P_r(2,1));
A(7,9) = -c*(R(1,1)/P_r(2,1));
A(7,11) = c*(R(3,3)*A1_r(2,1)/P_r(2,1)^2);
%partial derivatives with respect to y2 (right image)
A(8,1) = c*(R(2,2)/P_r(2,1));
A(8,2) = -c*(R(3,3)*A2_r(2,1)/P_r(2,1)^2);
A(8,3) = -(c/P_r(2,1)^2)*((R(1,1)*P_r(2,1)^2)+(R(1,1)*A2_r(2,1)^2));
A(8,4) = -c*(sin(kappa)*(1+(A2_r(2,1)^2)/P_r(2,1)^2)-cos(kappa)*A1_r(2,1)*A2_r(2,1)/P_r(2,1)^2);
A(8,5) = c*A1_r(2,1)/P_r(2,1);
A(8,10) = -c*R(2,2)/P_r(2,1);
A(8,11) = c*(R(3,3)*A2_r(2,1)/P_r(2,1)^2);
%partial derivatives with respect to x3 (left image)
A(9,12) = -c*(R(1,1)/P_l(3,1));
A(9,14) = c*(R(3,3)*A1_l(3,1)/P_l(3,1)^2);
%partial derivatives with respect to y3 (left image)
A(10,13) = -c*(R(2,2)/P_l(3,1));
A(10,14) = c*(R(3,3)*A2_l(3,1)/P_l(3,1)^2);
%partial derivatives with respect to x3 (right image)
A(11,2) = -c*(R(3,3)*A1_r(3,1)/P_r(3,1)^2);
A(11,3) = -c*(R(1,1)*A1_r(3,1)*A2_r(3,1)/P_r(3,1)^2);
A(11,4) = -c*((sin(kappa)*(A1_r(3,1)*A2_r(3,1)/P_r(3,1)^2))-cos(kappa)*(1+(A1_r(3,1)^2)/P_r(3,1)^2));
A(11,5) = -c*(A2_r(3,1)/P_r(3,1));
A(11,12) = -c*(R(1,1)/P_r(3,1));
A(11,14) = c*(R(3,3)*A1_r(3,1)/P_r(3,1)^2);
%partial derivatives with respect to y3 (right image)
A(12,1) = c*(R(2,2)/P_r(3,1));
A(12,2) = -c*(R(3,3)*A2_r(3,1)/P_r(3,1)^2);
A(12,3) = -(c/P_r(3,1)^2)*((R(1,1)*P_r(3,1)^2)+(R(1,1)*A2_r(3,1)^2));
A(12,4) = -c*(sin(kappa)*(1+(A2_r(3,1)^2)/P_r(3,1)^2)-cos(kappa)*A1_r(3,1)*A2_r(3,1)/P_r(3,1)^2);
A(12,5) = c*A1_r(3,1)/P_r(3,1);
A(12,13) = -c*R(2,2)/P_r(3,1);
A(12,14) = c*(R(3,3)*A2_r(3,1)/P_r(3,1)^2);
%partial derivatives with respect to x4 (left image)
A(13,15) = -c*(R(1,1)/P_l(4,1));
A(13,17) = c*(R(3,3)*A1_l(4,1)/P_l(4,1)^2);
%partial derivatives with respect to y4 (left image)
A(14,16) = -c*(R(2,2)/P_l(4,1));
A(14,17) = c*(R(3,3)*A2_l(4,1)/P_l(4,1)^2);
%partial derivatives with respect to x4 (right image)
A(15,2) = -c*(R(3,3)*A1_r(4,1)/P_r(4,1)^2);
A(15,3) = -c*(R(1,1)*A1_r(4,1)*A2_r(4,1)/P_r(4,1)^2);
A(15,4) = -c*((sin(kappa)*(A1_r(4,1)*A2_r(4,1)/P_r(4,1)^2))-cos(kappa)*(1+(A1_r(4,1)^2)/P_r(4,1)^2));
A(15,5) = -c*(A2_r(4,1)/P_r(4,1));
A(15,15) = -c*(R(1,1)/P_r(4,1));
A(15,17) = c*(R(3,3)*A1_r(4,1)/P_r(4,1)^2);
%partial derivatives with respect to y4 (right image)
A(16,1) = c*(R(2,2)/P_r(4,1));
A(16,2) = -c*(R(3,3)*A2_r(4,1)/P_r(4,1)^2);
A(16,3) = -(c/P_r(4,1)^2)*((R(1,1)*P_r(4,1)^2)+(R(1,1)*A2_r(4,1)^2));
A(16,4) = -c*(sin(kappa)*(1+(A2_r(4,1)^2)/P_r(4,1)^2)-cos(kappa)*A1_r(4,1)*A2_r(4,1)/P_r(4,1)^2);
A(16,5) = c*A1_r(4,1)/P_r(4,1);
A(16,16) = -c*R(2,2)/P_r(4,1);
A(16,17) = c*(R(3,3)*A2_r(4,1)/P_r(4,1)^2);
%partial derivatives with respect to x5 (left image)
A(17,18) = -c*(R(1,1)/P_l(5,1));
A(17,20) = c*(R(3,3)*A1_l(5,1)/P_l(5,1)^2);
%partial derivatives with respect to y5 (left image)
A(18,19) = -c*(R(2,2)/P_l(5,1));
A(18,20) = c*(R(3,3)*A2_l(5,1)/P_l(5,1)^2);
%partial derivatives with respect to x5 (right image)
A(19,2) = -c*(R(3,3)*A1_r(5,1)/P_r(5,1)^2);
A(19,3) = -c*(R(1,1)*A1_r(5,1)*A2_r(5,1)/P_r(5,1)^2);
A(19,4) = -c*((sin(kappa)*(A1_r(5,1)*A2_r(5,1)/P_r(5,1)^2))-cos(kappa)*(1+(A1_r(5,1)^2)/P_r(5,1)^2));
A(19,5) = -c*(A2_r(5,1)/P_r(5,1));
A(19,18) = -c*(R(1,1)/P_r(5,1));
A(19,20) = c*(R(3,3)*A1_r(5,1)/P_r(5,1)^2);
%partial derivatives with respect to y5 (right image)
A(20,1) = c*(R(2,2)/P_r(5,1));
A(20,2) = -c*(R(3,3)*A2_r(5,1)/P_r(5,1)^2);
A(20,3) = -(c/P_r(5,1)^2)*((R(1,1)*P_r(5,1)^2)+(R(1,1)*A2_r(5,1)^2));
A(20,4) = -c*(sin(kappa)*(1+(A2_r(5,1)^2)/P_r(5,1)^2)-cos(kappa)*A1_r(5,1)*A2_r(5,1)/P_r(5,1)^2);
A(20,5) = c*A1_r(5,1)/P_r(5,1);
A(20,19) = -c*R(2,2)/P_r(5,1);
A(20,20) = c*(R(3,3)*A2_r(5,1)/P_r(5,1)^2);
%partial derivatives with respect to x6 (left image)
A(21,21) = -c*(R(1,1)/P_l(6,1));
A(21,23) = c*(R(3,3)*A1_l(6,1)/P_l(6,1)^2);
%partial derivatives with respect to y6 (left image)
A(22,22) = -c*(R(2,2)/P_l(6,1));
A(22,23) = c*(R(3,3)*A2_l(6,1)/P_l(6,1)^2);
%partial derivatives with respect to x6 (right image)
A(23,2) = -c*(R(3,3)*A1_r(6,1)/P_r(6,1)^2);
A(23,3) = -c*(R(1,1)*A1_r(6,1)*A2_r(6,1)/P_r(6,1)^2);
A(23,4) = -c*((sin(kappa)*(A1_r(6,1)*A2_r(6,1)/P_r(6,1)^2))-cos(kappa)*(1+(A1_r(6,1)^2)/P_r(6,1)^2));
A(23,5) = -c*(A2_r(6,1)/P_r(6,1));
A(23,21) = -c*(R(1,1)/P_r(6,1));
A(23,23) = c*(R(3,3)*A1_r(6,1)/P_r(6,1)^2);
%partial derivatives with respect to y6 (right image)
A(24,1) = c*(R(2,2)/P_r(6,1));
A(24,2) = -c*(R(3,3)*A2_r(6,1)/P_r(6,1)^2);
A(24,3) = -(c/P_r(6,1)^2)*((R(1,1)*P_r(6,1)^2)+(R(1,1)*A2_r(6,1)^2));
A(24,4) = -c*(sin(kappa)*(1+(A2_r(6,1)^2)/P_r(6,1)^2)-cos(kappa)*A1_r(6,1)*A2_r(6,1)/P_r(6,1)^2);
A(24,5) = c*A1_r(6,1)/P_r(6,1);
A(24,22) = -c*R(2,2)/P_r(6,1);
A(24,23) = c*(R(3,3)*A2_r(6,1)/P_r(6,1)^2);

A(1:24,3:5) = A(1:24,3:5)/r;

AT = A';
N = AT*A;

N1 = inv(N);

%x from collinearity equation for the left image
Fx_l = zeros(6,1);
Fx_l(1,1) = -c*(A1_l(1,1)/P_l(1,1));
Fx_l(2,1) = -c*(A1_l(2,1)/P_l(2,1));
Fx_l(3,1) = -c*(A1_l(3,1)/P_l(3,1));
Fx_l(4,1) = -c*(A1_l(4,1)/P_l(4,1));
Fx_l(5,1) = -c*(A1_l(5,1)/P_l(5,1));
Fx_l(6,1) = -c*(A1_l(6,1)/P_l(6,1));
%y collinearity equation for the left image
Fy_l = zeros(6,1);
Fy_l(1,1) = -c*(A2_l(1,1)/P_l(1,1));
Fy_l(2,1) = -c*(A2_l(2,1)/P_l(2,1));
Fy_l(3,1) = -c*(A2_l(3,1)/P_l(3,1));
Fy_l(4,1) = -c*(A2_l(4,1)/P_l(4,1));
Fy_l(5,1) = -c*(A2_l(5,1)/P_l(5,1));
Fy_l(6,1) = -c*(A2_l(6,1)/P_l(6,1));
%x from collinearity equation for the right image
Fx_r = zeros(6,1);
Fx_r(1,1) = -c*(A1_r(1,1)/P_r(1,1));
Fx_r(2,1) = -c*(A1_r(2,1)/P_r(2,1));
Fx_r(3,1) = -c*(A1_r(3,1)/P_r(3,1));
Fx_r(4,1) = -c*(A1_r(4,1)/P_r(4,1));
Fx_r(5,1) = -c*(A1_r(5,1)/P_r(5,1));
Fx_r(6,1) = -c*(A1_r(6,1)/P_r(6,1));
%y from collinearity equation for the right image
Fy_r = zeros(6,1);
Fy_r(1,1) = -c*(A2_r(1,1)/P_r(1,1));
Fy_r(2,1) = -c*(A2_r(2,1)/P_r(2,1));
Fy_r(3,1) = -c*(A2_r(3,1)/P_r(3,1));
Fy_r(4,1) = -c*(A2_r(4,1)/P_r(4,1));
Fy_r(5,1) = -c*(A2_r(5,1)/P_r(5,1));
Fy_r(6,1) = -c*(A2_r(6,1)/P_r(6,1));
dl = zeros(24,1);
j = 1;
for i = 1:6
    dl(j,1) = left(i,1)-Fx_l(i,1);
    dl(j+1,1) = left(i,2)-Fy_l(i,1);
    dl(j+2,1) = right(i,1)-Fx_r(i,1);
    dl(j+3,1) = right(i,2)-Fy_r(i,1);
    j = j+4;
end
%matrix ×
dX = N1*AT*dl;
by = dX(1,1)+by;
bz = dX(2,1)+bz;
omega = (dX(3,1)/10000)+omega*200/pi();
phi = (dX(4,1)/10000)+phi*200/pi();
kappa = (dX(5,1)/10000)+kappa*200/pi();

U = A*dX - dl;

loop1=0;
for loop2 = 1:10
    if abs(dX(6,1))>0.01 || abs(dX(7,1))>0.01 || abs(dX(8,1))>0.01 || abs(dX(9,1))>0.01 || abs(dX(10,1))>0.01 || abs(dX(11,1))>0.01 || abs(dX(12,1))>0.01 || abs(dX(13,1))>0.01 || abs(dX(14,1))>0.01 || abs(dX(15,1))>0.01 || abs(dX(16,1))>0.01 || abs(dX(17,1))>0.01 || abs(dX(18,1))>0.01 || abs(dX(19,1))>0.01 || abs(dX(20,1))>0.01 || abs(dX(21,1))>0.01 || abs(dX(22,1))>0.01 || abs(dX(23,1))>0.01
        % Rùöê
        omega = omega*pi()/200;
        phi = phi*pi()/200;
        kappa = kappa*pi()/200;
        R = zeros(3,3);
        R(1,1) = cos(phi)*cos(kappa);
        R(1,2) = sin(omega)*sin(phi)*cos(kappa)+cos(omega)*sin(kappa);
        R(1,3) = -cos(omega)*sin(phi)*cos(kappa)+sin(omega)*sin(kappa);
        R(2,1) = -cos(phi)*sin(kappa);
        R(2,2) = -sin(omega)*sin(phi)*sin(kappa)+cos(omega)*cos(kappa);
        R(2,3) = cos(omega)*sin(phi)*sin(kappa)+sin(omega)*cos(kappa);
        R(3,1) = sin(phi);
        R(3,2) = -sin(omega)*cos(phi);
        R(3,3) = cos(omega)*cos(phi);
       
        f = zeros(3,3);
        f(1,1) = -sin(phi)*cos(kappa);
        f(1,2) = sin(omega)*cos(phi)*cos(kappa);
        f(1,3) = -cos(omega)*cos(phi)*cos(kappa);
        f(2,1) = sin(phi)*sin(kappa);
        f(2,2) = -sin(omega)*cos(phi)*sin(kappa);
        f(2,3) = cos(omega)*cos(phi)*sin(kappa);
        f(3,1) = cos(phi);
        f(3,2) = sin(omega)*sin(phi);
        f(3,3) = -cos(omega)*sin(phi);
        
        A1_l = zeros(6,1);
        A1_l(1,1) = R(1,1)*(geo(1,1)-X0l)+R(1,2)*(geo(1,2)-Y0l)+R(1,3)*(geo(1,3)-Z0l);
        A1_l(2,1) = R(1,1)*(geo(2,1)-X0l)+R(1,2)*(geo(2,2)-Y0l)+R(1,3)*(geo(2,3)-Z0l);
        A1_l(3,1) = R(1,1)*(geo(3,1)-X0l)+R(1,2)*(geo(3,2)-Y0l)+R(1,3)*(geo(3,3)-Z0l);
        A1_l(4,1) = R(1,1)*(geo(4,1)-X0l)+R(1,2)*(geo(4,2)-Y0l)+R(1,3)*(geo(4,3)-Z0l);
        A1_l(5,1) = R(1,1)*(geo(5,1)-X0l)+R(1,2)*(geo(5,2)-Y0l)+R(1,3)*(geo(5,3)-Z0l);
        A1_l(6,1) = R(1,1)*(geo(6,1)-X0l)+R(1,2)*(geo(6,2)-Y0l)+R(1,3)*(geo(6,3)-Z0l);
        
        A2_l = zeros(6,1);
        A2_l(1,1) = R(2,1)*(geo(1,1)-X0l)+R(2,2)*(geo(1,2)-Y0l)+R(2,3)*(geo(1,3)-Z0l);
        A2_l(2,1) = R(2,1)*(geo(2,1)-X0l)+R(2,2)*(geo(2,2)-Y0l)+R(2,3)*(geo(2,3)-Z0l);
        A2_l(3,1) = R(2,1)*(geo(3,1)-X0l)+R(2,2)*(geo(3,2)-Y0l)+R(2,3)*(geo(3,3)-Z0l);
        A2_l(4,1) = R(2,1)*(geo(4,1)-X0l)+R(2,2)*(geo(4,2)-Y0l)+R(2,3)*(geo(4,3)-Z0l);
        A2_l(5,1) = R(2,1)*(geo(5,1)-X0l)+R(2,2)*(geo(5,2)-Y0l)+R(2,3)*(geo(5,3)-Z0l);
        A2_l(6,1) = R(2,1)*(geo(6,1)-X0l)+R(2,2)*(geo(6,2)-Y0l)+R(2,3)*(geo(6,3)-Z0l);
       
        P_l = zeros(6,1);
        P_l(1,1) = R(3,1)*(geo(1,1)-X0l)+R(3,2)*(geo(1,2)-Y0l)+R(3,3)*(geo(1,3)-Z0l);
        P_l(2,1) = R(3,1)*(geo(2,1)-X0l)+R(3,2)*(geo(2,2)-Y0l)+R(3,3)*(geo(2,3)-Z0l);
        P_l(3,1) = R(3,1)*(geo(3,1)-X0l)+R(3,2)*(geo(3,2)-Y0l)+R(3,3)*(geo(3,3)-Z0l);
        P_l(4,1) = R(3,1)*(geo(4,1)-X0l)+R(3,2)*(geo(4,2)-Y0l)+R(3,3)*(geo(4,3)-Z0l);
        P_l(5,1) = R(3,1)*(geo(5,1)-X0l)+R(3,2)*(geo(5,2)-Y0l)+R(3,3)*(geo(5,3)-Z0l);
        P_l(6,1) = R(3,1)*(geo(6,1)-X0l)+R(3,2)*(geo(6,2)-Y0l)+R(3,3)*(geo(6,3)-Z0l);
       
        A1_r = zeros(6,1);
        A1_r(1,1) = R(1,1)*(geo(1,1)-bx)+R(1,2)*(geo(1,2)-by)+R(1,3)*(geo(1,3)-bz);
        A1_r(2,1) = R(1,1)*(geo(2,1)-bx)+R(1,2)*(geo(2,2)-by)+R(1,3)*(geo(2,3)-bz);
        A1_r(3,1) = R(1,1)*(geo(3,1)-bx)+R(1,2)*(geo(3,2)-by)+R(1,3)*(geo(3,3)-bz);
        A1_r(4,1) = R(1,1)*(geo(4,1)-bx)+R(1,2)*(geo(4,2)-by)+R(1,3)*(geo(4,3)-bz);
        A1_r(5,1) = R(1,1)*(geo(5,1)-bx)+R(1,2)*(geo(5,2)-by)+R(1,3)*(geo(5,3)-bz);
        A1_r(6,1) = R(1,1)*(geo(6,1)-bx)+R(1,2)*(geo(6,2)-by)+R(1,3)*(geo(6,3)-bz);
     
        A2_r = zeros(6,1);
        A2_r(1,1) = R(2,1)*(geo(1,1)-bx)+R(2,2)*(geo(1,2)-by)+R(2,3)*(geo(1,3)-bz);
        A2_r(2,1) = R(2,1)*(geo(2,1)-bx)+R(2,2)*(geo(2,2)-by)+R(2,3)*(geo(2,3)-bz);
        A2_r(3,1) = R(2,1)*(geo(3,1)-bx)+R(2,2)*(geo(3,2)-by)+R(2,3)*(geo(3,3)-bz);
        A2_r(4,1) = R(2,1)*(geo(4,1)-bx)+R(2,2)*(geo(4,2)-by)+R(2,3)*(geo(4,3)-bz);
        A2_r(5,1) = R(2,1)*(geo(5,1)-bx)+R(2,2)*(geo(5,2)-by)+R(2,3)*(geo(5,3)-bz);
        A2_r(6,1) = R(2,1)*(geo(6,1)-bx)+R(2,2)*(geo(6,2)-by)+R(2,3)*(geo(6,3)-bz);
   
        P_r = zeros(6,1);
        P_r(1,1) = R(3,1)*(geo(1,1)-bx)+R(3,2)*(geo(1,2)-by)+R(3,3)*(geo(1,3)-bz);
        P_r(2,1) = R(3,1)*(geo(2,1)-bx)+R(3,2)*(geo(2,2)-by)+R(3,3)*(geo(2,3)-bz);
        P_r(3,1) = R(3,1)*(geo(3,1)-bx)+R(3,2)*(geo(3,2)-by)+R(3,3)*(geo(3,3)-bz);
        P_r(4,1) = R(3,1)*(geo(4,1)-bx)+R(3,2)*(geo(4,2)-by)+R(3,3)*(geo(4,3)-bz);
        P_r(5,1) = R(3,1)*(geo(5,1)-bx)+R(3,2)*(geo(5,2)-by)+R(3,3)*(geo(5,3)-bz);
        P_r(6,1) = R(3,1)*(geo(6,1)-bx)+R(3,2)*(geo(6,2)-by)+R(3,3)*(geo(6,3)-bz);
       
        A = zeros(24,23);
       
        A(1,6) = -c*(R(1,1)/P_l(1,1));
        A(1,8) = c*(R(3,3)*A1_l(1,1)/P_l(1,1)^2);
       
        A(2,7) = -c*(R(2,2)/P_l(1,1));
        A(2,8) = c*(R(3,3)*A2_l(1,1)/P_l(1,1)^2);
       
        A(3,2) = -c*(R(3,3)*A1_r(1,1)/P_r(1,1)^2);
        A(3,3) = -c*(R(1,1)*A1_r(1,1)*A2_r(1,1)/P_r(1,1)^2);
        A(3,4) = -c*((sin(kappa)*(A1_r(1,1)*A2_r(1,1)/P_r(1,1)^2))-cos(kappa)*(1+(A1_r(1,1)^2)/P_r(1,1)^2));
        A(3,5) = -c*(A2_r(1,1)/P_r(1,1));
        A(3,6) = -c*(R(1,1)/P_r(1,1));
        A(3,8) = c*(R(3,3)*A1_r(1,1)/P_r(1,1)^2);
       
        A(4,1) = c*(R(2,2)/P_r(1,1));
        A(4,2) = -c*(R(3,3)*A2_r(1,1)/P_r(1,1)^2);
        A(4,3) = -(c/P_r(1,1)^2)*((R(1,1)*P_r(1,1)^2)+(R(1,1)*A2_r(1,1)^2));
        A(4,4) = -c*(sin(kappa)*(1+(A2_r(1,1)^2)/P_r(1,1)^2)-cos(kappa)*A1_r(1,1)*A2_r(1,1)/P_r(1,1)^2);
        A(4,5) = c*A1_r(1,1)/P_r(1,1);
        A(4,7) = -c*R(2,2)/P_r(1,1);
        A(4,8) = c*(R(3,3)*A2_r(1,1)/P_r(1,1)^2);
       
        A(5,9) = -c*(R(1,1)/P_l(2,1));
        A(5,11) = c*(R(3,3)*A1_l(2,1)/P_l(2,1)^2);
       
        A(6,10) = -c*(R(2,2)/P_l(2,1));
        A(6,11) = c*(R(3,3)*A2_l(2,1)/P_l(2,1)^2);
       
        A(7,2) = -c*(R(3,3)*A1_r(2,1)/P_r(2,1)^2);
        A(7,3) = -c*(R(1,1)*A1_r(2,1)*A2_r(2,1)/P_r(2,1)^2);
        A(7,4) = -c*((sin(kappa)*(A1_r(2,1)*A2_r(2,1)/P_r(2,1)^2))-cos(kappa)*(1+(A1_r(2,1)^2)/P_r(2,1)^2));
        A(7,5) = -c*(A2_r(2,1)/P_r(2,1));
        A(7,9) = -c*(R(1,1)/P_r(2,1));
        A(7,11) = c*(R(3,3)*A1_r(2,1)/P_r(2,1)^2);
      
        A(8,1) = c*(R(2,2)/P_r(2,1));
        A(8,2) = -c*(R(3,3)*A2_r(2,1)/P_r(2,1)^2);
        A(8,3) = -(c/P_r(2,1)^2)*((R(1,1)*P_r(2,1)^2)+(R(1,1)*A2_r(2,1)^2));
        A(8,4) = -c*(sin(kappa)*(1+(A2_r(2,1)^2)/P_r(2,1)^2)-cos(kappa)*A1_r(2,1)*A2_r(2,1)/P_r(2,1)^2);
        A(8,5) = c*A1_r(2,1)/P_r(2,1);
        A(8,10) = -c*R(2,2)/P_r(2,1);
        A(8,11) = c*(R(3,3)*A2_r(2,1)/P_r(2,1)^2);
       
        A(9,12) = -c*(R(1,1)/P_l(3,1));
        A(9,14) = c*(R(3,3)*A1_l(3,1)/P_l(3,1)^2);
        
        A(10,13) = -c*(R(2,2)/P_l(3,1));
        A(10,14) = c*(R(3,3)*A2_l(3,1)/P_l(3,1)^2);
       
        A(11,2) = -c*(R(3,3)*A1_r(3,1)/P_r(3,1)^2);
        A(11,3) = -c*(R(1,1)*A1_r(3,1)*A2_r(3,1)/P_r(3,1)^2);
        A(11,4) = -c*((sin(kappa)*(A1_r(3,1)*A2_r(3,1)/P_r(3,1)^2))-cos(kappa)*(1+(A1_r(3,1)^2)/P_r(3,1)^2));
        A(11,5) = -c*(A2_r(3,1)/P_r(3,1));
        A(11,12) = -c*(R(1,1)/P_r(3,1));
        A(11,14) = c*(R(3,3)*A1_r(3,1)/P_r(3,1)^2);
       
        A(12,1) = c*(R(2,2)/P_r(3,1));
        A(12,2) = -c*(R(3,3)*A2_r(3,1)/P_r(3,1)^2);
        A(12,3) = -(c/P_r(3,1)^2)*((R(1,1)*P_r(3,1)^2)+(R(1,1)*A2_r(3,1)^2));
        A(12,4) = -c*(sin(kappa)*(1+(A2_r(3,1)^2)/P_r(3,1)^2)-cos(kappa)*A1_r(3,1)*A2_r(3,1)/P_r(3,1)^2);
        A(12,5) = c*A1_r(3,1)/P_r(3,1);
        A(12,13) = -c*R(2,2)/P_r(3,1);
        A(12,14) = c*(R(3,3)*A2_r(3,1)/P_r(3,1)^2);
        
        A(13,15) = -c*(R(1,1)/P_l(4,1));
        A(13,17) = c*(R(3,3)*A1_l(4,1)/P_l(4,1)^2);
        
        A(14,16) = -c*(R(2,2)/P_l(4,1));
        A(14,17) = c*(R(3,3)*A2_l(4,1)/P_l(4,1)^2);
       
        A(15,2) = -c*(R(3,3)*A1_r(4,1)/P_r(4,1)^2);
        A(15,3) = -c*(R(1,1)*A1_r(4,1)*A2_r(4,1)/P_r(4,1)^2);
        A(15,4) = -c*((sin(kappa)*(A1_r(4,1)*A2_r(4,1)/P_r(4,1)^2))-cos(kappa)*(1+(A1_r(4,1)^2)/P_r(4,1)^2));
        A(15,5) = -c*(A2_r(4,1)/P_r(4,1));
        A(15,15) = -c*(R(1,1)/P_r(4,1));
        A(15,17) = c*(R(3,3)*A1_r(4,1)/P_r(4,1)^2);
        
        A(16,1) = c*(R(2,2)/P_r(4,1));
        A(16,2) = -c*(R(3,3)*A2_r(4,1)/P_r(4,1)^2);
        A(16,3) = -(c/P_r(4,1)^2)*((R(1,1)*P_r(4,1)^2)+(R(1,1)*A2_r(4,1)^2));
        A(16,4) = -c*(sin(kappa)*(1+(A2_r(4,1)^2)/P_r(4,1)^2)-cos(kappa)*A1_r(4,1)*A2_r(4,1)/P_r(4,1)^2);
        A(16,5) = c*A1_r(4,1)/P_r(4,1);
        A(16,16) = -c*R(2,2)/P_r(4,1);
        A(16,17) = c*(R(3,3)*A2_r(4,1)/P_r(4,1)^2);
       
        A(17,18) = -c*(R(1,1)/P_l(5,1));
        A(17,20) = c*(R(3,3)*A1_l(5,1)/P_l(5,1)^2);
       
        A(18,19) = -c*(R(2,2)/P_l(5,1));
        A(18,20) = c*(R(3,3)*A2_l(5,1)/P_l(5,1)^2);
        
        A(19,2) = -c*(R(3,3)*A1_r(5,1)/P_r(5,1)^2);
        A(19,3) = -c*(R(1,1)*A1_r(5,1)*A2_r(5,1)/P_r(5,1)^2);
        A(19,4) = -c*((sin(kappa)*(A1_r(5,1)*A2_r(5,1)/P_r(5,1)^2))-cos(kappa)*(1+(A1_r(5,1)^2)/P_r(5,1)^2));
        A(19,5) = -c*(A2_r(5,1)/P_r(5,1));
        A(19,18) = -c*(R(1,1)/P_r(5,1));
        A(19,20) = c*(R(3,3)*A1_r(5,1)/P_r(5,1)^2);
        
        A(20,1) = c*(R(2,2)/P_r(5,1));
        A(20,2) = -c*(R(3,3)*A2_r(5,1)/P_r(5,1)^2);
        A(20,3) = -(c/P_r(5,1)^2)*((R(1,1)*P_r(5,1)^2)+(R(1,1)*A2_r(5,1)^2));
        A(20,4) = -c*(sin(kappa)*(1+(A2_r(5,1)^2)/P_r(5,1)^2)-cos(kappa)*A1_r(5,1)*A2_r(5,1)/P_r(5,1)^2);
        A(20,5) = c*A1_r(5,1)/P_r(5,1);
        A(20,19) = -c*R(2,2)/P_r(5,1);
        A(20,20) = c*(R(3,3)*A2_r(5,1)/P_r(5,1)^2);
        
        A(21,21) = -c*(R(1,1)/P_l(6,1));
        A(21,23) = c*(R(3,3)*A1_l(6,1)/P_l(6,1)^2);
        
        A(22,22) = -c*(R(2,2)/P_l(6,1));
        A(22,23) = c*(R(3,3)*A2_l(6,1)/P_l(6,1)^2);
        
        A(23,2) = -c*(R(3,3)*A1_r(6,1)/P_r(6,1)^2);
        A(23,3) = -c*(R(1,1)*A1_r(6,1)*A2_r(6,1)/P_r(6,1)^2);
        A(23,4) = -c*((sin(kappa)*(A1_r(6,1)*A2_r(6,1)/P_r(6,1)^2))-cos(kappa)*(1+(A1_r(6,1)^2)/P_r(6,1)^2));
        A(23,5) = -c*(A2_r(6,1)/P_r(6,1));
        A(23,21) = -c*(R(1,1)/P_r(6,1));
        A(23,23) = c*(R(3,3)*A1_r(6,1)/P_r(6,1)^2);
        
        A(24,1) = c*(R(2,2)/P_r(6,1));
        A(24,2) = -c*(R(3,3)*A2_r(6,1)/P_r(6,1)^2);
        A(24,3) = -(c/P_r(6,1)^2)*((R(1,1)*P_r(6,1)^2)+(R(1,1)*A2_r(6,1)^2));
        A(24,4) = -c*(sin(kappa)*(1+(A2_r(6,1)^2)/P_r(6,1)^2)-cos(kappa)*A1_r(6,1)*A2_r(6,1)/P_r(6,1)^2);
        A(24,5) = c*A1_r(6,1)/P_r(6,1);
        A(24,22) = -c*R(2,2)/P_r(6,1);
        A(24,23) = c*(R(3,3)*A2_r(6,1)/P_r(6,1)^2);

        A(1:24,3:5) = A(1:24,3:5)/r;
        
        AT = A';
        N = AT*A;
        
        N1 = inv(N);
        
        Fx_l = zeros(6,1);
        Fx_l(1,1) = -c*(A1_l(1,1)/P_l(1,1));
        Fx_l(2,1) = -c*(A1_l(2,1)/P_l(2,1));
        Fx_l(3,1) = -c*(A1_l(3,1)/P_l(3,1));
        Fx_l(4,1) = -c*(A1_l(4,1)/P_l(4,1));
        Fx_l(5,1) = -c*(A1_l(5,1)/P_l(5,1));
        Fx_l(6,1) = -c*(A1_l(6,1)/P_l(6,1));
       
        Fy_l = zeros(6,1);
        Fy_l(1,1) = -c*(A2_l(1,1)/P_l(1,1));
        Fy_l(2,1) = -c*(A2_l(2,1)/P_l(2,1));
        Fy_l(3,1) = -c*(A2_l(3,1)/P_l(3,1));
        Fy_l(4,1) = -c*(A2_l(4,1)/P_l(4,1));
        Fy_l(5,1) = -c*(A2_l(5,1)/P_l(5,1));
        Fy_l(6,1) = -c*(A2_l(6,1)/P_l(6,1));
        
        Fx_r = zeros(6,1);
        Fx_r(1,1) = -c*(A1_r(1,1)/P_r(1,1));
        Fx_r(2,1) = -c*(A1_r(2,1)/P_r(2,1));
        Fx_r(3,1) = -c*(A1_r(3,1)/P_r(3,1));
        Fx_r(4,1) = -c*(A1_r(4,1)/P_r(4,1));
        Fx_r(5,1) = -c*(A1_r(5,1)/P_r(5,1));
        Fx_r(6,1) = -c*(A1_r(6,1)/P_r(6,1));
        
        Fy_r = zeros(6,1);
        Fy_r(1,1) = -c*(A2_r(1,1)/P_r(1,1));
        Fy_r(2,1) = -c*(A2_r(2,1)/P_r(2,1));
        Fy_r(3,1) = -c*(A2_r(3,1)/P_r(3,1));
        Fy_r(4,1) = -c*(A2_r(4,1)/P_r(4,1));
        Fy_r(5,1) = -c*(A2_r(5,1)/P_r(5,1));
        Fy_r(6,1) = -c*(A2_r(6,1)/P_r(6,1));
        dl = zeros(24,1);
        j = 1;
        for i = 1:6
            dl(j,1) = left(i,1)-Fx_l(i,1);
            dl(j+1,1) = left(i,2)-Fy_l(i,1);
            dl(j+2,1) = right(i,1)-Fx_r(i,1);
            dl(j+3,1) = right(i,2)-Fy_r(i,1);
            j = j+4;
        end
  
        dX = N1*AT*dl;
        by = dX(1,1)+by;
        bz = dX(2,1)+bz;
        omega = (dX(3,1)/10000)+omega*200/pi();
        phi = (dX(4,1)/10000)+phi*200/pi();
        kappa = (dX(5,1)/10000)+kappa*200/pi();
        
        U = A*dX - dl;
        loop1 = loop1 + 1;
    end
end
UT = U';
s02 = (UT*U)/1;

Vx = s02*N1;

sby = sqrt(Vx(1,1));
sbz = sqrt(Vx(2,2));
somega = sqrt(Vx(3,3))/10000;
sphi = sqrt(Vx(4,4))/10000;
skappa = sqrt(Vx(5,5))/10000;

sX1 = sqrt(Vx(6,6));
sY1 = sqrt(Vx(7,7));
sZ1 = sqrt(Vx(8,8));
sX2 = sqrt(Vx(9,9));
sY2 = sqrt(Vx(10,10));
sZ2 = sqrt(Vx(11,11));
sX3 = sqrt(Vx(12,12));
sY3 = sqrt(Vx(13,13));
sZ3 = sqrt(Vx(14,14));
sX4 = sqrt(Vx(15,15));
sY4 = sqrt(Vx(16,16));
sZ4 = sqrt(Vx(17,17));
sX5 = sqrt(Vx(18,18));
sY5 = sqrt(Vx(19,19)); 
sZ5 = sqrt(Vx(20,20));
sX6 = sqrt(Vx(21,21));
sY6 = sqrt(Vx(22,22));
sZ6 = sqrt(Vx(23,23));