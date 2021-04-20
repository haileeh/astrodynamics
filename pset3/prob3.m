%% prob 3 inefficient interplanetary transfers
% scale up transfer orbit dV from energy minimizing HTO
% consider range of angles when scaled TO dV = 1.5*V_h

%G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
%mEarth = 5.974*10^24; %[kg]
%rEarth = 6380; %[km] Earth radius
%rMars = 3400; %[km] Mars radius
%muEarth = 5.16*10^12; %[km^3/hr^2]
%muEarth_s = 3.987*10^5; %[km^3/s^2]
%muMars = 5.57*10^11; %[km^3/hr^2]
%muMars_s = 4.298*10^4; %[km^3/s^2]
%soi_Earth = 0.929*10^6; %[km]
%soi_Mars = 0.578*10^6;%[km]
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
rEarth_Sun = 149.60*10^6; %[km]
rMars_Sun = 227.9*10^6; %[km]

% values that may be used
v_Earth = sqrt(muSun_s/rEarth_Sun); %[km/sec]
aTrans = (rEarth_Sun + rMars_Sun)/2; %[km]
v_Mars = sqrt(muSun_s/rMars_Sun); %km/sec
%v_transE = sqrt(muSun_s * (2/rEarth_Sun - 1/aTrans));  %this is v_pi
%% Earth to Mars
% Values in km and sec
v_hohmann = sqrt(muSun_s/rEarth_Sun)*(sqrt(2*rMars_Sun/(rEarth_Sun+rMars_Sun))-1); %vpi - v_earth 
dv = 1.5*v_hohmann;
p.mu = muSun_s;

% tspan
tf = 2*pi*aTrans^(3/2) / sqrt(p.mu);
tspan = linspace(0,tf,1000);
opts = odeset;

theta = linspace(-54,54,20); %deg, with respect to y
x = zeros(length(tspan),length(theta)); y = x;
for i=1:length(theta)
    inits = [rEarth_Sun; 0; dv*sind(theta(i)); v_Earth + dv*cosd(theta(i))];
    [t,z] = ode45(@rhs,tspan,inits,opts,p); 
    x(:,i) = z(:,1);
    y(:,i) = z(:,2);
end
% circle of Earth's orbit about Sun
th = 0:pi/50:2*pi;
xunit = rEarth_Sun * cos(th);
yunit = rEarth_Sun * sin(th);
% circle of Mars's orbit about Sun
xunitM = rMars_Sun * cos(th);
yunitM = rMars_Sun * sin(th);

figure;
plot(0,0,'yo','LineWidth',2); hold on;
plot(xunit,yunit,'b','LineWidth',2); 
plot(xunitM,yunitM,'r','LineWidth',2);
plot(x,y);%,'--g','LineWidth',2);
axis equal
title('Earth to Mars Trajectories');
xlabel('X Position [km]');
ylabel('Y Position [km]');

%% Mars to Earth
% Values in km and sec
v_hohmann = sqrt(muSun_s/rMars_Sun)*(1-sqrt(2*rEarth_Sun/(rEarth_Sun+rMars_Sun))); %v_mars - v_pi 
dv = 1.5*v_hohmann;
p.mu = muSun_s;

% tspan
tf = 2*pi*aTrans^(3/2) / sqrt(p.mu);
tspan = linspace(0,tf,1000);
opts = odeset;

theta = linspace(-51,51,20); %deg, with respect to y
x = zeros(length(tspan),length(theta)); y = x;
for i=1:length(theta)
    inits = [rMars_Sun; 0; -dv*sind(theta(i)); v_Mars - dv*cosd(theta(i))];
    [t,z] = ode45(@rhs,tspan,inits,opts,p); 
    x(:,i) = z(:,1);
    y(:,i) = z(:,2);
end


figure;
plot(0,0,'yo','LineWidth',2); hold on;
plot(xunit,yunit,'b','LineWidth',2); 
plot(xunitM,yunitM,'r','LineWidth',2);
plot(x,y);%,'--g','LineWidth',2);
axis equal
title('Mars to Earth Trajectories')
xlabel('X Position [km]');
ylabel('Y Position [km]');

keyboard
%% non-dimensional Earth to Mars
% find normalized mu
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
muSun_hr = muSun_s * 3600^2; %[km^3/hr^2]
muSun_day = muSun_hr * 24^2; %[km^3/day^2]
mu_n = muSun_day / rEarth_Sun^3; %1/days^2
r2 = rMars_Sun/rEarth_Sun;
v_hohmann_n = sqrt(mu_n)*(sqrt(2*r2 /(1+r2))-1); %[1/day]

% solve for orbit
p.mu = mu_n; % mu is in 1/days^2
theta = -54; % 54 to -54
dv_n = 1.5*v_hohmann_n; % normalized value
vE_n = v_Earth/rEarth_Sun*(60*60*24);

inits = [rEarth_Sun/rEarth_Sun; 0; dv_n*sind(theta); vE_n + dv_n*cosd(theta)];

opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
a = aTrans / rEarth_Sun;
tf = 2*pi*a^(3/2) / sqrt(p.mu);
tspan = linspace(0,tf,100);

[t,z] = ode45(@rhs,tspan,inits,opts,p);
x = z(:,1);
y = z(:,2);

% circle of Earth's orbit about Sun
th = 0:pi/50:2*pi;
xunit = rEarth_Sun/rEarth_Sun * cos(th);
yunit = rEarth_Sun/rEarth_Sun * sin(th);
% circle of Mars's orbit about Sun
xunitM = rMars_Sun/rEarth_Sun * cos(th);
yunitM = rMars_Sun/rEarth_Sun * sin(th);

figure;
plot(0,0,'yo','LineWidth',2); hold on;
plot(xunit,yunit,'b','LineWidth',2); 
plot(xunitM,yunitM,'r','LineWidth',2);
plot(x,y,'--g','LineWidth',2);
axis equal

%% Mars to Earth non-dimensional
r2 = rMars_Sun/rEarth_Sun;
v_hohmann_n = sqrt(mu_n/r2)*(1-sqrt(2*1 /(1+r2))); %[1/day]
%v_hohmann = sqrt(muSun_s/rMars_Sun)*(1-sqrt(2*rEarth_Sun/(rEarth_Sun+rMars_Sun))); %v_mars - v_pi 


% solve for orbit
p.mu = mu_n; % mu is in 1/days^2
theta = -51; % 51 to -51
dv_n = 1.5*v_hohmann_n; % normalized value
vM_n = v_Mars/rEarth_Sun*(60*60*24);

inits = [rMars_Sun/rEarth_Sun; 0; -dv_n*sind(theta); vM_n - dv_n*cosd(theta)];

opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
a = aTrans / rEarth_Sun;
tf = 2*pi*a^(3/2) / sqrt(p.mu);
tspan = linspace(0,tf,100);

[t,z] = ode45(@rhs,tspan,inits,opts,p);
x = z(:,1);
y = z(:,2);

figure;
plot(0,0,'yo','LineWidth',2); hold on;
plot(xunit,yunit,'b','LineWidth',2); 
plot(xunitM,yunitM,'r','LineWidth',2);
plot(x,y,'--g','LineWidth',2);
axis equal