function [x,y,t] = FandGFcns(R0,V0,mu)

% norms
V_norm = norm(V0);
R_norm =  norm(R0);

% eccentricity vector
e_vec = cross(V0,cross(R0,V0))/mu - R0/R_norm;
e = norm(e_vec);

% find semi-major axis
a = (2/R_norm - V_norm^2/mu)^-1;

% find angular momentum
h_vec = cross(R0,V0);
h = norm(h_vec);

% set up coordinate system
ihatE = e_vec/e;
ihatH = h_vec/h;
ihatP = cross(ihatE,ihatH);

A = [ihatE, ihatP, ihatH]; 
r0 = A*R0;
v0 = A*V0;

% find E0
r0_norm = norm(r0);
v0_norm = norm(v0);
sigma0 = dot(r0,v0) / sqrt(mu);
E0 = atan2((sigma0/sqrt(a)),(1-r0_norm/a)); %equivalent expression
%E0 = acos((a-r0_norm)/(a*e));

% solve for F,G
n = sqrt(mu/a^3);
t0 = (1/n) * (E0 - e*sin(E0));
Et = linspace(E0,E0+(2*pi),100)';
Etilde = Et-E0;
t = (Et - e*sin(Et))/n;
F = 1 - (a/r0_norm)*(1-cos(Etilde));
G = (t-t0) + sqrt(a^3/mu)*(sin(Etilde)-Etilde);
%Fdot = -sqrt(mu*a)/(r*r0_norm)*sin(Etilde);
%Gdot = 1 - (a/r)*(1-cos(Etilde));
x0 = r0(1);
y0 = r0(2);
xdot0 = v0(1);
ydot0 = v0(2);
x = x0*F + xdot0*G;
y = y0*F + ydot0*G;

% figure; 
% plot(x,y,'--k','LineWidth',3); hold on;
% xlabel('X-Position');
% ylabel('Y-Position');
% %xlim([-0.2,1.8]); ylim([-1,1]);
% title('Orbit from F,G Functions')

% figure
% subplot(2,1,1);
% plot(x);
% subplot(2,1,2);
% plot(y)
