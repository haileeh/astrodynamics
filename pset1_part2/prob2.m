%% prob 2.1 - Molniya Orbit
% orbital period
T = 0.5*(23*3600 + 56 * 60 + 4.091); %[sec]
inc = 63.4 * pi/180; %[rads]
rEarth = 6371; %[km]
rAlpha = 1.1*rEarth; %[km]
G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
mEarth = 5.974*10^24; %[kg]
mu =  G*mEarth; %[ km^3 / s^2 ]
a = (T*sqrt(mu)/(2*pi))^(2/3); %[ km ]
mage = 1 - rAlpha/a; 
param = a*(1 - mage^2);
vAlpha = (1 + mage)*sqrt(mu/param); %[km/s]

% boundary angles of E
Eb1 = acos(mage); 
Eb2 = 2*pi - Eb1;
% boundary angles of M
Mb1 = Eb1 - mage*sin(Eb1);
Mb2 = Eb2 - mage*sin(Eb2);

% find percentage of time spent above the equator
nMotion = sqrt(mu/a^3);
t1 = Mb1/nMotion;
t2 = Mb2/nMotion;
deltaT = t2-t1;
propT = deltaT/T; %[s]

%% 2.2 - plot the orbit in Cartesian and Polar
r0 = [rAlpha; 0; 0];
v0 = [0; vAlpha; 0];

% plot the orbit
p.mu_AUD = mu;
inits = [rAlpha 0 0 vAlpha];
% time span (one period)
tspan = linspace(0,T+0.05*T,300);
% solve ode
opts = odeset;
[tstar,z]=ode45(@rhs,tspan,inits,opts,p);
x = z(:,1);
y = z(:,2);

figure
plot(z(:,1),z(:,2)); hold on;
title('X-Position vs. Y-Position from Cartesian');
xlim([-5e4,1e4]); ylim([-2e4,2e4]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', 'k', 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', 'k', 'LineStyle', linestyle);

% time histories for r and theta
r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);

figure;
subplot(2,1,1);
plot(tstar,r); hold on; grid on;
title('Time history of r from Cartesian');
xlabel('Time'); ylabel('r(t)');
subplot(2,1,2);
plot(tstar,theta);hold on; grid on;
title('Time history of \theta from Cartesian');
xlabel('Time'); ylabel('\theta(t)');

%%% repeat using polar coordinates
r0_polar = norm(r0);
theta0 = atan2(r0(2),r0(1));
rdot0 = 0;
thetadot0 = vAlpha/rAlpha;
inits = [r0_polar theta0 rdot0 thetadot0];

[tz,z]=ode45(@rhs_polar,tspan,inits,opts,p);
r = z(:,1);
theta = z(:,2);
x_polar = r.*cos(theta);
y_polar = r.*sin(theta);

figure
plot(x_polar,y_polar); hold on;
title('X-Position vs. Y-Position from Polar');
xlim([-5e4,1e4]); ylim([-2e4,2e4]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', 'k', 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', 'k', 'LineStyle', linestyle);

figure;
subplot(2,1,1);
plot(tz,r); hold on; grid on;
title('Time history of r from Polar');
xlabel('Time'); ylabel('r(t)');
subplot(2,1,2);
plot(tz,theta);hold on; grid on;
title('Time history of \theta from Polar');
xlabel('Time'); ylabel('\theta(t)');

%% cross equator time
idx = find(theta <= pi/2 + 0.04);
idx = find(theta(idx) >= pi/2 - 0.04);
tz(idx)
percDiff = (tz(idx)-t1)/t1
