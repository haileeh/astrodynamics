% richardson's analytical solution - he is going off of L1
% not sure why ydot0 is so large

% constants
n_13 = 3; % not sure what this is
mEarth = 5.9722*10^24; %[kg]
mSun = 1.989*10^30; %[kg]
%mu = mEarth/(mSun+mEarth);
mu = 3.040357143*10^-6; % from Richardson; mu = G*E in this context; not sure what E is
au = 149.5978714*10^6; %km %R12
LU = au;
TU = 5.030853327*10^3; % seconds; from Richardson's
n = 1.990986606*10^-7; % rad/s
G = (6.67408e-11); % m^3 / kg / s^2
ne = -sqrt((mSun+mEarth)*G/((au)*10^3)^3); % rad/s
Az = 110000; %km % from Richardson's paper? nominal orbit 

% could be more precise below
rL = 1.497610042*10^6; %rE in richardson %km % distance of secondary body from the Lagrange point (1)
Az = Az/rL;
gammaL = rL /au;
c2 = (1/gammaL^3)* (mu + (1-mu)*gammaL^3/(1-gammaL)^3);
c3 = (1/gammaL^3)* (mu - (1-mu)*gammaL^4/(1-gammaL)^4);
c4 = (1/gammaL^3)* (mu + (1-mu)*gammaL^5/(1-gammaL)^5);
% linearized frequency lambda by solving
% syms lambdas real
% eqn = lambdas^4 + (c2-2)*lambdas^2 - (c2-1)*(1+2*c2);
% sol = solve(eqn, lambdas);
% if sol(1) > 0
%     lambda = sol(1);
% else
%     lambda = sol(2);
% end

lambda2 = (-(c2-2) + sqrt((c2-2)^2 + 4*1*(c2-1)*(1+2*c2)))/(2*1);
lambda = sqrt(lambda2);
k = 2*lambda / (lambda^2 + 1 - c2);

d1 = 3*lambda^2 / k * (k*(6*lambda^2 -1) - 2*lambda);
d2 = 8*lambda^2 / k * (k*(11*lambda^2 -1) - 2*lambda);

a21 = 3*c3*(k^2-2)/ (4*(1+2*c2));
a22 = 3*c3 / (4*(1+2*c2));
a23 = -3*c3*lambda/(4*k*d1) * (3*k^3*lambda - 6*k*(k-lambda)+4);
a24 = -3*c3*lambda/(4*k*d1) * (2+3*k*lambda);

d21 = -c3/(2*lambda^2);
d31 = 3/(64*lambda^2)*(4*c3*a24+c4);
d32 = 3/(64*lambda^2) * (4*c3*(a23-d21) + c4*(4+k^2));

b21 = -3*c3*lambda/(2*d1)*(3*k*lambda-4);
b22 = 3*c3*lambda / d1;
b31 = 3/(8*d2)* (8*lambda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2)) + (9*lambda^2+1+2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2)));
b32 = 1/d2* (9*lambda*(c3*(k*b22 + d21 - 2*a24) - c4) + (3/8)*(9*lambda^2 + 1 + 2*c2)*(4*c3*(k*a24-b22)+k*c4));

a31 = -9*lambda/(4*d2) * (4*c3*(k*a23-b21) + k*c4*(4+k^2)) + (9*lambda^2+1-c2)/(2*d2)*(3*c3*(2*a23-k*b21)+c4*(2+3*k^2));
a32 = -1/d2*( 9*lambda/4 * (4*c3*(k*a24-b22) + k*c4) + 3/2*(9*lambda^2+1-c2)*(c3*(k*b22+d21-2*a24)-c4));
% frequency corrections s1 and s2
s1 = 1/(2*lambda*(lambda*(1+k^2)-2*k)) * (3/2*c3*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21) - 3/8*c4*(3*k^4-8*k^2+8));
s2 = 1/(2*lambda*(lambda*(1+k^2)-2*k)) * (3/2*c3*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+5*d21) + 3/8*c4*(12-k^2));

a1 = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-k^2);
a2 = 3/2*c3*(a24-2*a22) + 9/8*c4;

l1 = a1+2*lambda^2*s1;
l2 = a2+2*lambda^2*s2;
Delta = lambda^2-c2; % nomalized units
Delta_unNorm = Delta * LU^2;

%Ax = sqrt((-l2*Az^2-Delta_unNorm)/l1); % works with unNorm
Ax = -206000 / rL; 
Ay = k*Ax;
% s1_unnorm = s1*LU^2;
% s2_unnorm = s2*LU^2;
% omega = 1 + s1_unnorm*Ax^2 + s2_unnorm*Az^2
omega = 1 + s1*Ax^2 + s2*Az^2; % of the correct order O(A_x^n)
% solve
t = 0; % don't use syms
s = n*t;
phi = 0; %determines the family of orbits
tau = omega*s;
tau1 = lambda*tau + phi; %n=[rad/s], t =s, omega = [km^2]?

x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Ax^2)*cos(2*tau1) + ...
    (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);
y = k*Ax*sin(tau1) + (b21*Ax^2-b22*Az^2)*sin(2*tau1)+ (b31*Ax^3-b32*Ax*Az^2)*sin(3*tau1);
deln = 2 - n_13;
z = deln*(Az*cos(tau1) + d21*Ax*Az*(cos(2*tau1) - 3) + (d32*Az*Ax^2 -d31*Az^3)*cos(3*tau1));

% initial conditions for a halo orbit in an L1 centered coordinate frame
% normalized by the Lagrange point secondary distance (rL)
tstar = vpasolve(y,t);
x0 = subs(x,t,tstar);
y0 = subs(y,t,tstar);
z0 = subs(z,t,tstar); %wrong
ydot0 = k*Ax*(lambda*cos(tau1)) + (b21*Ax^2 - b22*Az^2)*(2*lambda*cos(2*tau1)) + ...
    (b31*Ax^3 - b32*Ax*Az^2)*(3*lambda*cos(3*tau1));
ydot0 = subs(ydot0,t,tstar); % wrong

% convert to CR3BP frame
% x0_gammaL = x0*gammaL;
% z0_gammaL = z0*gammaL;
% ydot0_gammaL = ydot0*gammaL;
% 
% x0 = (LU - rL+x0_gammaL)/LU;
% z0 = z0_gammaL/LU;
% ydot0 = ydot0_gammaL*TU/LU;

T = 2*pi/(lambda*omega*n); % normalized units
%omega_Az = s1*((-l2*Az^2-Delta_unNorm)/l1) + s2*Az^2;