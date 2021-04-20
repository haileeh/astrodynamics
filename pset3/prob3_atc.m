%% prob 3 inefficient interplanetary transfers
% scale up transfer orbit dV from energy minimizing HTO
% consider range of angles when scaled TO dV = 1.5*V_h
G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
mEarth = 5.974*10^24; %[kg]
rEarth = 6380; %[km] Earth radius
rMars = 3400; %[km] Mars radius
muEarth = 5.16*10^12; %[km^3/hr^2]
muEarth_s = 3.987*10^5; %[km^3/s^2]
muMars = 5.57*10^11; %[km^3/hr^2]
muMars_s = 4.298*10^4; %[km^3/s^2]
soi_Earth = 0.929*10^6; %[km]
soi_Mars = 0.578*10^6;%[km]
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
rEarth_Sun = 149.60*10^6; %[km]
rMars_Sun = 227.9*10^6; %[km]

% values that may be used
v_Earth = sqrt(muSun_s/rEarth_Sun); %[km/sec]
aTrans = (rEarth_Sun + rMars_Sun)/2; %[km]
v_Mars = sqrt(muSun_s/rMars_Sun); %km/sec
v_transE = sqrt(muSun_s * (2/rEarth_Sun - 1/aTrans));  %this is v_pi
%% Earth to Mars
% Values in km and sec
v_hohmann = sqrt(muSun_s/rEarth_Sun)*(sqrt(2*rMars_Sun/(rEarth_Sun+rMars_Sun))-1);
dv = 1.5*v_hohmann;
p.mu = muSun_s;

% tspan
tf = 2*pi*aTrans^(3/2) / sqrt(p.mu);
tspan = linspace(0,tf,1000);
opts = odeset;
% inits
theta = 54; %deg % now vary this; %theta defined with respect to y, 54
inits = [rEarth_Sun; 0; dv*sind(theta); v_Earth+dv*cosd(theta)];
[t,z] = ode45(@rhs,tspan,inits,opts,p); 

x = z(:,1);
y = z(:,2);

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
plot(x,y,'--g','LineWidth',2);
axis equal

function zdot = rhs(t,z,p)

xdot = -p.mu*z(1)/((z(1)^2+z(2)^2))^(3/2);
ydot = -p.mu*z(2)/((z(1)^2+z(2)^2))^(3/2);
zdot = [z(3);z(4);xdot;ydot];

end