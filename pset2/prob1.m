%% prob 1.1 Hyperbolic transfer orbit - escape from Earth into Hohmann
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

r1p = 200; %[km]

v_eo = sqrt(muEarth_s*(1/r1p)); % answer
%v_soi = sqrt(muEarth_s*1/soi_Earth); %v2

v_pi = sqrt(muEarth_s*(2/r1p - 2/(r1p+soi_Earth)));
v1inf = v_pi-v_eo;

% energy is constant equation
v1p = sqrt(v1inf^2 + 2*muEarth_s/r1p);

deltaV_esc1 = v1p - v_eo;  % final answer

%% prob 1.2 hyperbolic capture by Mars from HTO; then parking orbit
% velocity should be smaller going from inner to outer planet
v_alpha = sqrt(muEarth_s*(2/soi_Earth - 2/(r1p+soi_Earth)));

r2p = 200; %[km]
v_mo = sqrt(muMars_s*(1/r2p)); % answer
v2inf = v_alpha-v_mo; 
% energy is constant equation
v2p = sqrt(v2inf^2 + 2*muMars_s/r2p); %answer
deltaV_cap1 = v2p - v_mo;

% what happens if v2p is zero?
r_far = 2*muMars_s/v2inf^2; %[km]
v_far = sqrt(muMars_s/r_far);
r = 2*muMars_s/(v_mo^2-v2inf^2)
r_inf = muMars_s/v2inf^2
% none of the above are correct
%% 1.3 hyperbolic escape from Mars to HTO
r3p = 200; %[km]
v_mo =sqrt(muMars_s*(1/r3p)); %answer
v_pi = sqrt(muMars_s*(2/r3p - 2/(r3p+soi_Mars)));
v3inf = v_pi-v_mo;

% energy is constant equation
v3p = sqrt(v3inf^2 + 2*muMars_s/r3p); %answer
deltaV_esc2 = v3p - v_mo;  % final answer
%% 1.4 - hyperbolic capture by Earth from HTO
v_alpha = sqrt(muMars_s*(2/soi_Mars - 2/(r3p+soi_Mars)));
r4p = 200;

v_eo = sqrt(muEarth_s*(1/r4p)); % answer
v4inf = v_alpha-v_eo 
% energy is constant equation
v4p = sqrt(v4inf^2 + 2*muEarth_s/r4p) %answer

deltaV_cap2 = v4p - v_eo;

%% 1.5 home sweet home
deltaV_tot = deltaV_esc1 + deltaV_cap1 + deltaV_esc2 + deltaV_cap2;

%rocket equation
c = 4500; % exhaust velocity

deltaM_dryMass = exp(deltaV_tot/c)-1; %convert to percent

% now solve for Delta for each trajectory
% first solve for beta
%beta = acos(1/(1+rp*v_inf^2/mu)); %mu is planet you are leaving
% then solve for length of otherside ?
d1 = r1p*v1p/v1inf %[km]
d2 = r2p*v2p/v2inf %[km]
d3 = r3p*v3p/v3inf %[km]
d4 = r4p*v4p/v4inf %[km]

% dist Earth to Sun = 149.6 million km