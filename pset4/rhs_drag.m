function zdot = rhs_drag(t,z,p)

% xdot = -p.mu_AUD*z(1)/((z(1)^2+z(2)^2))^(3/2);
% ydot = -p.mu_AUD*z(2)/((z(1)^2+z(2)^2))^(3/2);
radius = z;
alt = radius/1000 - p.R_J; %km
rho = p.rho0 * exp(-alt/p.H); %kg / m^3
dadt = -sqrt(p.mu_Jupiter*radius)*rho*p.effAe; % m/s

zdot = dadt;