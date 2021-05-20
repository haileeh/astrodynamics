% prob 1 - SSO
rho_deg_yr = 360;
rho_deg_sec = rho_deg_yr / (365*24*60*60);
rho = rho_deg_sec * pi / 180; % rad/s, precession rate
%rho = 1.99096871*10^-7; %rad/s, precession rate
mu = 398600.440; % km3/s2
R_E = 6378; %km
J2 = 1.08263*10^-3; %second zonal term
alt = 3300; %km
a = alt + R_E; %km
e = 0;
inc1 = acosd(-2*rho*a^(7/2)*(1-e^2)^2/(3*J2*R_E^2*sqrt(mu)));
%%  part 2 how does eccentricity change affect argument of perigee
dwdt = -(1.5 * (sqrt(mu)*J2*R_E^2/( (1-e^2)^2*a^(7/2)) ) ) * (5/2*sind(inc1)^2 - 2); % rad/sec
dwdt1 = dwdt * 180/pi * (24*60*60) % deg/day

% now, eccenctiricty e = 0.25, how does this change argument of perigee?
e = 0.25;
inc2 = acosd(-2*rho*a^(7/2)*(1-e^2)^2/(3*J2*R_E^2*sqrt(mu)))
dwdt = -(1.5 * (sqrt(mu)*J2*R_E^2/( (1-e^2)^2*a^(7/2)) ) ) * (5/2*sind(inc2)^2 - 2); % rad/sec
dwdt2 = dwdt * 180/pi * (24*60*60) % deg/day

% plot argument of perigee, assuming everything else fixed
T = 2*pi*sqrt(a^3 / mu); %seconds
T_min = T / 60; 
T_hr = T_min / 60;
T_days = T_hr / 24;
tspan = linspace(0,10*T_days,100);
w1 = dwdt1*tspan; % deg
w2 = dwdt2*tspan; % deg

figure;
plot(tspan, w1, 'r', 'LineWidth',2);
hold on; grid on;
plot(tspan, w2, '--g', 'LineWidth',2);
xlabel('Time [days]');
ylabel('\omega [deg]');
xlim([0 tspan(end)]);
