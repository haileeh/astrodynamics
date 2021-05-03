function [L1,L2,L3] = getLpoints(mu)
% return normalized value
% mEarth = 5.9722*10^24; %[kg]
% mSun = 1.989*10^30; %[kg]
% au = 149.5978714*10^6; %km
%
% M1 = mSun; M2 = mEarth; R = au;
%
% L1 = (M2/(3*M1))^(1/3);
% L2 = L1; %on other side of smaller body
% L3 = (7/12)*M2/M1;

syms func x y z
R1 = sqrt((x+mu)^2 + y^2 + z^2);
R2 = sqrt((x+mu-1)^2 + y^2 + z^2);
func = (1/2)*(x^2+y^2) + (1-mu)/R1 + mu/R2;
dx = diff(func, x);
yl = 0; zl = 0;
dx_l = subs(dx,   y, yl);
dx_l = subs(dx_l, z, zl);

L1 = double( vpasolve( dx_l == 0,  0) );
L2 = double( vpasolve( dx_l == 0,  1) );
L3 = double( vpasolve( dx_l == 0, -1) );
