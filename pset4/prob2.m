% prob 2
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
dT = 1; % time step in seconds
T = 0; % start time in seconds
it = 0; % iteration counter

%% solve for non-changing values
radius = (R_J+alt)*1000; %m
effAe = CD*A/m_sc;
% find period
P = 2*pi*sqrt(radius^3/mu_Jupiter); % sec
Pm = P / 60; % mins
mm = 1440 / Pm; % mean motion
%% loop
% while it < 1e4
%     % atmospheric density
%     rho = rho0*exp(-alt/H); %kg/m^3
%     % decrement in Orbital period
%     dT_sec = dT;%*3600;
%     dP = 3*pi*radius*rho*(A*CD/m_sc)*dT_sec; % sec
%     
%     %decay = dP / P * mm / dT; % rev/day/day
%     
%     P = P - dP;
%     T = T + dT;
%     % update orbital radius
%     radius = (mu_Jupiter * P^2/ (4*pi^2) )^(1/3);
%     alt = radius/1000 - R_J; % km
%     % solve for new alt
%     if alt <= alt_f % km
%         fprintf('Time to decay: %.2f sec \n',T);
%         break
%     end
%     it = it+1;
% end
% keyboard
%% loop 2
dT = 0.01; % min
dT_sec = dT*60;
while it < 1e4
    % atmospheric density
    rho = rho0*exp(-alt/H); %kg/m^3
    % decrement in Orbital period
    
    %da = -sqrt(mu_Jupiter*radius)*rho*effAe* dT_sec; % m
    da = -sqrt(mu_Jupiter)*rho*effAe * dT_sec / 2; % sqrt([m])

    % update orbital radius
    %radius = radius + da;
    radius = (sqrt(radius) + da)^2;
    T = T + dT;

    alt = radius/1000 - R_J; % km
    % solve for new alt
    if alt <= alt_f % km
        fprintf('Time to decay: %.2f min \n',T);
        break
    end
    it = it+1;
end
%% loop 3
% dT = 1; % min
% dT_sec = dT*60;
% Ph = Pm / 60;
% Pd = Ph / 24;
% while it < 1e4
%     % atmospheric density
%     rho = rho0*exp(-alt/H); %kg/m^3
%     % decrement in Orbital period
%     
%     %da = -sqrt(mu_Jupiter*radius)*rho*effAe* dT_sec; % m
%     deltaA_rev = -2*pi*effAe*rho*radius^2;
%     %da = -sqrt(mu_Jupiter)*rho*effAe * dT_sec / 2; % sqrt([m])
% 
%     % update orbital radius
%     radius = radius + deltaA_rev;
%     %radius = (sqrt(radius) + da)^2;
%     T = T + Pd;
% 
%     alt = radius/1000 - R_J; % km
%     % solve for new alt
%     if alt <= alt_f % km
%         fprintf('Time to decay: %.2f day \n',T);
%         break
%     end
%     it = it+1;
% end