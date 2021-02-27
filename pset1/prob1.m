% prob 1
r_Earth = 6.37812*10^3; %[km]
m_Earth = 5.974*10^24; %[kg]
r_Earth_Earth = 0;

r_Moon = 1.738*10^3; %[km]
m_Moon = 7.3483*10^22; %[kg]
r_Moon_Earth = 3.844*10^5; %[km]

G = 6.67259*10^-11; % [ m^3 kg^-1 s^-2]

mu = G*(m_Earth + m_Moon);

r_CM_Earth = m_Moon*r_Moon_Earth / (m_Moon + m_Earth);

r_Moon_CM = r_Moon_Earth - r_CM_Earth;
r_Moon_CM_m = r_Moon_CM*10^3;

% period of motion
T_sec = 2*pi*r_Moon_CM_m^(3/2) / sqrt(mu);
T_day = T_sec/(60*60*24);

% velocities
v_Moon = 2*pi*r_Moon_CM_m / T_sec; %[m/s]
v_Earth = 2*pi*r_CM_Earth*10^3 / T_sec; % [m/s]

%% prob 3
% use km and sec
%mu_km = mu * (1/1000)^3;
G_km = G * (1/1000)^3;
mu_Earth = G_km * m_Moon;
mu_Moon = G_km * m_Earth;

rEarthInit = r_CM_Earth;%-4670.8 ; % km % positive
rMoonInit = r_Moon_CM;%3.7973*10^5; % km
vMoonInit = 1.0309; % km/s
vEarthInit = -0.01268; % km/s

% inits = [rE0, rM0, vE0, vM0]
inits = [rEarthInit, 0, rMoonInit, 0, 0, vEarthInit, 0, vMoonInit];

p.mu_Earth = mu_Earth;
p.mu_Moon = mu_Moon;

% time span
%tf = 24*3600*30;%24*3600*365;
tf = T_sec;
tspan = linspace(0,tf,300);

% solve ode
opts = odeset;
[t,z]=ode45(@rhs_2body,tspan,inits,opts,p);

figure
plot(z(:,3)-z(:,1),z(:,4)-z(:,2));
hold on;
title('Combined Motion: X-Position vs. Y-Position');
xlim([-5e5, 5e5]);
ylim([-5e5, 5e5]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);

%% prob 3b
tOrbit = T_sec/(2*pi);
rMoonInitND = rMoonInit / r_Moon_Earth;
rEarthInitND = rEarthInit / r_Moon_Earth;
vMoonInitND = vMoonInit/r_Moon_Earth * tOrbit;
vEarthInitND = vEarthInit / r_Moon_Earth * tOrbit;
T = 2*pi;
G = G_km*(tOrbit^2)/ (r_Moon_Earth^3); %1/kg
%muEND = m_Earth*G;
%muMND = m_Moon*G; 
% non-dimensionalized
p.mu_Earth = m_Moon*G;
p.mu_Moon = m_Earth*G;

inits = [rEarthInitND, 0, rMoonInitND, 0, 0, vEarthInitND, 0, vMoonInitND];

tspan = linspace(0,T,300);

% solve ode
[t,z]=ode45(@rhs_2body,tspan,inits,opts,p);

figure
plot(z(:,3)-z(:,1),z(:,4)-z(:,2));
hold on;
title('Nondimensionalized Combined Motion: X-Position vs. Y-Position');
% xlim([-5e5, 5e5]);
% ylim([-5e5, 5e5]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);