%% prob2 v2
% modeling atmospheric drag
alt = 500; %km
R_J = 71492; %km - equatorial radius
m_sc = 150; % kg
A = 2; %m^2 surface area
CD = 4; % coeff of drag
Temp = 165; % K
H = 27; %km, scale height
molarMass = 2.016; % g/mol
Rstar = 8.31446261815324; % m^3 Pa/ (K mol)
p0 = 10^5; % Pa 
rho0 = p0*molarMass / (Rstar * Temp); % g / m^3
rho0 = rho0 / 1000; %kg/m^3; density at 1 bar (https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html)

alt_f = 400; %km
mu_Jupiter = 1.26686534*10^17; %m^3/s^2

radius = (R_J+alt)*1000; %m
effAe = CD*A/m_sc;

%% drag force
%nu = sqrt(mu / radius);
%aD = -effAe * rho * nu^2 / 2

% alt = radius/1000 - R_J; %km
% rho = rho0 * exp(-alt/H); %kg / m^3
% dadt = -sqrt(mu_Jupiter*radius)*rho*effAe; % m/s
inits = radius;
p.H = H;
p.R_J = R_J;
p.rho0 = rho0;
p.effAe = effAe;
p.mu_Jupiter = mu_Jupiter;
tspan = linspace(0,150,50);
opts = odeset;
[t,z] = ode45(@rhs_drag,tspan,inits,opts,p); % z in m
z_km = z / 1000; 
alt_t = z_km - R_J;
figure;
plot(t,alt_t);
% same result as other method
