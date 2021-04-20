% Hohmann transfer Earth to Mars
% Constants
%G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
%mEarth = 5.974*10^24; %[kg]
%rEarth = 6380; %[km]
%rMars = 3400; %[km]
% muEarth = 5.16*10^12; %[km^3/hr^2]
% muEarth_s = 3.987*10^5; %[km^3/s^2]
% muMars = 5.57*10^11; %[km^3/hr^2]
% muMars_s = 4.298*10^4; %[km^3/s^2]
%soi_Earth = 0.929*10^6; %[km]
%soi_Mars = 0.578*10^6;%[km]
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
muSun_hr = muSun_s * 3600^2; %[km^3/hr^2]
muSun_day = muSun_hr * 24^2; %[km^3/day^2]
rEarth_Sun = 149.60*10^6; %[km]
rMars_Sun = 227.9*10^6; %[km]
%alt = 200; % [km]
% new constants
%muSun_day2 = muSun_day / rEarth_Sun^3; % this matches sqrtMuSun
mu = 1; % how to get this from muSun_s
sqrtMuSun = sqrt(2.96*10^-4);

% find transfer time in days
%tH = pi*sqrt((r1+r2)^3/(8*mu));
aH = (rEarth_Sun+rMars_Sun)/2; % [km]
%tH = 0.5*sqrt(4*pi^2*aH^3/muSun_day); %[days]

% this in non-dimensional
aH_nd = aH/rEarth_Sun;
tH_nd = 0.5*sqrt(4*pi^2*aH_nd^3/mu);
tH_days = tH_nd / sqrtMuSun; %[days]

%% circular orbits - at what angle relative to Earth does Mars have to be located at launch
% such that a rendezvous will occur at the end of the Hohmann transfer?
r2 = rMars_Sun/rEarth_Sun; % non-dimensional
omega2 = sqrt(mu/r2^3); % ND
phi_E = pi; 
phi_M = omega2*tH_nd;
% Eqn 8.7 Curtis
alpha = phi_E - phi_M; %[rads]
alpha_deg = alpha * 180/pi; %[deg]

%% Provide three possible future launch dates. What will be the relevant 
% arrival dates for each transfer? Create a plot of launch date vs arrival 
% date for transfers between 2021 and 2031, from both Earth to Mars and Mars to Earth.
tH_yr = tH_days / 365.25;
% last transfer day
yr_LeavePrev = 2020.4679;
yr_ArrPrev = yr_LeavePrev + tH_yr;

% What is the synodic period? Mars orbit relative to Earth
r1 = rEarth_Sun/rEarth_Sun;
omega1 = sqrt(mu/r1^3); % ND
T_E = 2*pi/omega1; % more precisely, T_E = 365.26 days
T_M = 2*pi/omega2;% more precisely, T_M = 687.99 days
T_syn = T_E*T_M / (abs(T_E-T_M)); % ND
T_syn_days = T_syn / sqrtMuSun;
T_syn_years = T_syn_days / 365.25;

% Earth to Mars!
leaveDates = yr_LeavePrev + T_syn_years*(0:1:10);
yrLeave = floor(leaveDates);
dayLeave = 365.25*(leaveDates-yrLeave);
arriveDates = leaveDates + tH_yr;
yrArr = floor(arriveDates);
dayArr = 365.25*(arriveDates-yrArr);

fprintf('First possible launch date: Year: %d, Day: %.3f \n',yrLeave(2),dayLeave(2));
fprintf('First possible arrival date: Year: %d, Day: %.3f \n',yrArr(2),dayArr(2));
fprintf('Second possible launch date: Year: %d, Day: %.3f \n',yrLeave(3),dayLeave(3));
fprintf('Second possible arrival date: Year: %d, Day: %.3f \n',yrArr(3),dayArr(3));
fprintf('Third possible launch date: Year: %d, Day: %.3f \n',yrLeave(4),dayLeave(4));
fprintf('Third possible arrival date: Year: %d, Day: %.3f \n \n',yrArr(4),dayArr(4));

%% find Mars to Earth launch dates and arrival dates
% Curtis Example 8.2
t12 = tH_days;
nEarth = 2*pi/(T_E/sqrtMuSun); %[rad/day]
nMars = 2*pi/(T_M/sqrtMuSun); %[rad/day]

% phase angle between earth and Mars when spacecraft reaches Mars is 
phi_f = pi - nEarth*t12; %[rad]
N = 1;
t_wait = ( -2*phi_f - 2*pi*N ) / (nMars - nEarth); % [days]
leaveMars = arriveDates + t_wait/365.25;
yrL = floor(leaveMars);
dayL = 365.25*(leaveMars-yrL);
arriveEarth = leaveMars + tH_yr; 
yrA = floor(arriveEarth);
dayA =  365.25*(arriveEarth-yrA);
% fprintf('1. Leave Mars on year: %d, day: %.3f \n',yrL(1),dayL(1));
% fprintf('1. Arrive Earth on year: %d, day: %.3f \n',yrA(1),dayA(1));
% fprintf('2. Leave Mars on year: %d, day: %.3f \n',yrL(2),dayL(2));
% fprintf('2. Arrive Earth on year: %d, day: %.3f \n',yrA(2),dayA(2));
% fprintf('3. Leave Mars on year: %d, day: %.3f \n',yrL(3),dayL(3));
% fprintf('3. Arrive Earth on year: %d, day: %.3f \n',yrA(3),dayA(3));

% plot launch date vs arrival date
figure;
plot(leaveDates,arriveDates,'x','LineWidth',2); hold on; grid on;
plot(leaveMars,arriveEarth,'o','LineWidth',2);
xlabel('Departure Date [yr]');
ylabel('Arrival Date [yr]');
legend('Earth to Mars','Mars to Earth','Location','southeast');
xlim([2021 2031]);
ylim([2021 2031]);
title('Earth-Mars Launch and Arrival Schedule')