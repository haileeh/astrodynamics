% CR3BP
%G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
G = 1; %nondimensionalized
mEarth = 5.974*10^24; %[kg]
mMoon = 7.348*10^22; %[kg]
M1 = mEarth; M2 = mMoon;
m = M1+M2;
L = 384400; %[km]
OMEGA = sqrt(G*(M1+M2)/L^3);
% find positions relative to CM
% this puts CM at zero
mRatioE = M1/(M1+M2);
% r1 = -(1-mRatioE)*L;
% r2 = mRatioE*L;

mu = M2/(M1+M2);
mu1 = 1-mu; mu2 = mu;
% in the inertial frame,
% X1 = -mu2*cos(t); Y1 = -mu2*sin(t); Z1 = 0; % larger mass
% X2 = mu1*cos(t); Y2 = mu1*sin(t); Z2 = 0;
x1 = -mu2; x2 = mu1; % rotating frame

% r1 = sqrt((X+mu2*cos(t))^2+(Y+mu2*sin(t))^2+Z^2);
% r2 = sqrt((X-mu1*cos(t))^2+(Y-mu1*sin(t))^2+Z^2);
% r1 and r2 in rotating coordinates:
% syms x y z mu2x mu1x real
% r1 = sqrt((x+mu2x)^2+(y+mu2x)^2+z^2);
% r2 = sqrt((x-mu1x)^2+(y-mu1x)^2+z^2);
% % define gravitational potential; r1, r2 are distances of P from m1,m2
% U = -mu1x/r1 - mu2x/r2 - 0.5*mu1x*mu2x; 
% Ubar = -0.5*(x^2+y^2)+U;
% Ubar_x = simplify(diff(Ubar,x));
% Ubar_y = simplify(diff(Ubar,y));
% Ubar_z = simplify(diff(Ubar,z));
% Lagrangian in rotating frame is time indep
%L_rot = 0.5*((xdot-y)^2+(ydot+x)^2 + zdot^2)-U;

% call ode
inits = [1;0;0;0;0.01;0];
tspan = linspace(0,2*pi,50);
opts = odeset;
p.mu1 = mu1; p.mu2 = mu2;
[t,zed]=ode45(@cr3bp_ode,tspan,inits,opts,p);
x = zed(:,1);y = zed(:,2);z = zed(:,3);
vx = zed(:,4); vy = zed(:,5); vz = zed(:,6);
r1 = sqrt((x+mu2).^2+(y+mu2).^2+z.^2);
r2 = sqrt((x-mu1).^2+(y-mu1).^2+z.^2);
% define gravitational potential; r1, r2 are distances of P from m1,m2
U = -mu1/r1 - mu2/r2 - 0.5*mu1*mu2; 
for i=1:length(t)
    Ubar(i) = -0.5*(x(i)^2+y(i)^2) + U(i);
% energy integral of motion
    E(i) = 0.5*(vx(i)^2 + vy(i)^2 + vz(i)^2) + Ubar(i);
end



function zeddot = cr3bp_ode(t,zed,p)
x = zed(1); y = zed(2); z = zed(3);
xdot = zed(4); ydot = zed(5); zdot = zed(6);
mu1 = p.mu1;
mu2 = p.mu2;

Ubar_x = (mu1*(mu2 + x))/((mu2 + x)^2 + (mu2 + y)^2 + z^2)^1.5000 - x - (mu2*(mu1 - x))/((mu1 - x)^2 + (mu1 - y)^2 + z^2)^1.5000;
Ubar_y = (mu1*(mu2 + y))/((mu2 + x)^2 + (mu2 + y)^2 + z^2)^1.5000 - y - (mu2*(mu1 - y))/((mu1 - x)^2 + (mu1 - y)^2 + z^2)^1.5000;
Ubar_z = mu2*z/((mu1-x)^2+(mu1-y)^2+z^2)^1.5 + mu1*z/((mu2+x)^2+(mu2+y)^2+z^2)^1.5;
xddot = 2*ydot - Ubar_x;
yddot = -2*xdot - Ubar_y;
zddot = -Ubar_z;

zeddot = [xdot;ydot;zdot;xddot;yddot;zddot];
end