%% prob 3.1 - Molynia orbit
omega_deg = 270; % argument of periapsis
omega_rad = omega_deg * pi/180;
inc_deg = 63.4; % inclination
inc_rad = inc_deg*pi/180;
Omega_rad = pi/2; % longitude of ascending node
T_sec = 0.5*(23*3600 + 56 * 60 + 4.091); %[sec]
rEarth = 6371; %[km]
rAlpha = 1.1*rEarth; %[km]

G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
mEarth = 5.974*10^24; %[kg]
mu_km =  G*mEarth; %[ km^3 / s^2 ]
a = (T_sec*sqrt(mu_km)/(2*pi))^(2/3); %[ km ]
mage = 1 - rAlpha/a; % eccentricity
param = a*(1 - mage^2);
vAlpha = (1 + mage)*sqrt(mu_km/param); %[km/s]

r = [rAlpha; 0; 0];
v = [0; vAlpha; 0];

p.mu_AUD = mu_km;
inits = [rAlpha 0 0 vAlpha];
tspan = linspace(0,2*T_sec,1000);
% solve ode
opts = odeset;
[t,z]=ode45(@rhs,tspan,inits,opts,p);
X3 = z(:,1);
Y3 = z(:,2);
%[X3,Y3,t] = FandGFcns([rAlpha;0;0],[0;vAlpha;0],mu_km); % fix special case
groundTracks(Omega_rad,inc_rad,omega_rad,X3,Y3,t);

%% prob 3.2 geostationary orbit - really a geosynchronous orbit since nonzero inclination
omega_deg = 0; % argument of periapsis 
omega_rad = omega_deg * pi/180;
inc_deg = 10; % inclination
inc_rad = 10 * pi/180;
Omega_rad = 0; % longitude of ascending node
T_sec = 1436*60; % period of motion
mage = 0.1; % eccentricity
a = 42164; %[km]
rAlpha = a*(1 - mage); %[km]
param = a*(1 - mage^2);
vAlpha = (1 + mage)*sqrt(mu_km/param); %[km/s]

inits = [rAlpha 0 0 vAlpha];
tspan = linspace(0,T_sec+10,100);
% solve ode
opts = odeset;
[t,z]=ode45(@rhs,tspan,inits,opts,p);
X3 = z(:,1);
Y3 = z(:,2);

groundTracks(Omega_rad,inc_rad,omega_rad,X3,Y3,t);
