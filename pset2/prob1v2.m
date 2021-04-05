% Constants from problem
G = (6.67259*10^-11)/(1000^3); %[km^3 / kg s^2 ]
mEarth = 5.974*10^24; %[kg]
rEarth = 6380; %[km]
rMars = 3400; %[km]
muEarth = 5.16*10^12; %[km^3/hr^2]
muEarth_s = 3.987*10^5; %[km^3/s^2]
muMars = 5.57*10^11; %[km^3/hr^2]
muMars_s = 4.298*10^4; %[km^3/s^2]
soi_Earth = 0.929*10^6; %[km]
soi_Mars = 0.578*10^6;%[km]
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
rEarth_Sun = 149.60*10^6; %[km]
rMars_Sun = 227.9*10^6; %[km]
c = 4500; % [m/s]? exhaust velocity
c = c/1000; %[km/s]
alt = 200; % [km]

% values used in multiple sections
v_Earth = sqrt(muSun_s/rEarth_Sun); %[km/sec]
aTrans = (rEarth_Sun + rMars_Sun)/2; %[km]
v_Mars = sqrt(muSun_s/rMars_Sun); %km/sec
%% prob 1.1 - Hyperbolic transfer orbit - escape from Earth into Hohmann
r1p = rEarth+alt; %[km]
v_eo = sqrt(muEarth_s/r1p); %[km\sec]

v_transE = sqrt(muSun_s * (2/rEarth_Sun - 1/aTrans));  %this is v_pi

vinf_E = v_transE-v_Earth;  

% from energy equation
v1p = sqrt(vinf_E^2 + 2*muEarth_s/r1p);
deltaV_1 = abs(v1p - v_eo);

% print out answers
fprintf('Problem 1.1 \n');
fprintf('Orbital velocity, v_eo, is %.3f km/s \n',v_eo);
fprintf('Velocity of Transfer at Earth (Periapsis) is %.3f km/s\n',v_transE);
fprintf('V @ infinity, v_inf, is %.3f km/s\n',vinf_E);
fprintf('Velocity at Periapsis, v_p, is %.3f km/s\n',v1p);
fprintf('Delta V = %.3f km/s\n', deltaV_1);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');
%% 1.2 hyperbolic capture by Mars
% velocity should be smaller going from inner to outer planet
v_transM = sqrt(muSun_s*(2/rMars_Sun - 1/aTrans));   % this is v_alpha

r2p = rMars+alt; %[km]
v_mo = sqrt(muMars_s/r2p); % answer

vinf_M = v_Mars - v_transM;

% energy is constant equation
v2p = sqrt(vinf_M^2 + 2*muMars_s/r2p);

deltaV_2 = abs(v2p - v_mo);

fprintf('Problem 1.2 \n');
fprintf('Orbital velocity, v_mo, is %.3f km/s \n',v_mo);
fprintf('Velocity of Transfer at Mars (Periapsis) is %.3f km/s\n',v_transM);
fprintf('V @ infinity, v_inf, is %.3f km/s\n',vinf_M);
fprintf('Velocity at Periapsis, v_p, is %.3f km/s\n',v2p);
fprintf('Delta V = %.3f km/s\n', deltaV_2);

% zero deltav? v_mo = v2p, flyby
ecc = r2p*v2p^2/muMars_s - 1;
delta = 2*asin(1/ecc); % turning angle

fprintf('No delta-V, e = %.3f, delta = %.3f \n',ecc,delta);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

%% 1.3 escape Mars to HTO
r3p = r2p; %[km]

% same as v_transM
v_transM2 = sqrt(muSun_s * (2/rMars_Sun - 1/aTrans));  %this is v_pi
vinf_M2 = v_transM2-v_Mars; 

% energy equation
v3p = sqrt(vinf_M2^2 + 2*muMars_s/r3p);  

deltaV_3 = abs(v3p - v_mo);

% print out answers
fprintf('Problem 1.3 \n');
fprintf('Orbital velocity, v_mo, is %.3f km/s \n',v_mo);
fprintf('Velocity of Transfer at Earth (Periapsis) is %.3f km/s\n',v_transM2);
fprintf('V @ infinity, v_inf, is %.3f km/s\n',vinf_M2);
fprintf('Velocity at Periapsis, v_p, is %.3f km/s\n',v3p);
fprintf('Delta V = %.3f km/s\n', deltaV_3);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

%% 1.4 capture by Earth
% velocity should be larger going from outer to inter planet
v_transE2 = sqrt(muSun_s*(2/rEarth_Sun - 1/aTrans));   % this is v_alpha

r4p = r1p; %[km]

vinf_E2 = v_Earth-v_transE2;

% energy is constant equation
v4p = sqrt(vinf_E2^2 + 2*muEarth_s/r4p);

deltaV_4= abs(v4p - v_eo);

% print out answers
fprintf('Problem 1.4 \n');
fprintf('Orbital velocity, v_eo, is %.3f km/s \n',v_eo);
fprintf('Velocity of Transfer at Earth (Periapsis) is %.3f km/s\n',v_transE2);
fprintf('V @ infinity, v_inf, is %.3f km/s\n',vinf_E2);
fprintf('Velocity at Periapsis, v_p, is %.3f km/s\n',v4p);
fprintf('Delta V = %.3f km/s\n', deltaV_4);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

%% 1.5 home sweet home
deltaV_tot = deltaV_1 + deltaV_2 + deltaV_3 + deltaV_4;

deltaM_dryMass = (1-exp(-deltaV_tot/c))*100; %convert to percent

% now solve for Delta for each trajectory
% eqn 8.57
d1 = r1p*sqrt(1+2*muEarth_s/(r1p*vinf_E^2)); %[km]
d2 = r2p*sqrt(1+2*muMars_s/(r2p*vinf_M^2)); %[km]
d3 = r3p*sqrt(1+2*muMars_s/(r3p*vinf_M2^2)); %[km]
d4 = r4p*sqrt(1+2*muEarth_s/(r4p*vinf_E2^2)); %[km]

% compare: astronomical units
d1_AU = d1/rEarth_Sun;
d2_AU = d2/rEarth_Sun;
d3_AU = d3/rEarth_Sun;
d4_AU = d4/rEarth_Sun;