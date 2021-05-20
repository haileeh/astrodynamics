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
rho0 = rho0 / 1000; %kg/m^3; 
alt_f = 400; %km
mu_Jupiter = 1.26686534*10^17; %m^3/s^2
T = 0; % start time in seconds
it = 1; % iteration counter
dT_sec = 0.01;
radius = (R_J+alt)*1000; %m
effAe = CD*A/m_sc;

%% loop
while it < 1e5
    % atmospheric density
    rho(it) = rho0*exp(-alt(it)/H); %kg/m^3
    
    % decrement in semi-major axis
    da = -sqrt(mu_Jupiter*radius)*rho(it)*effAe* dT_sec; % m
    %da = -sqrt(mu_Jupiter)*rho*effAe * dT_sec / 2; % sqrt([m])
    
    it = it+1;
    % update orbital radius
    radius = radius + da;
    %radius = (sqrt(radius) + da)^2;
    T = T + dT_sec;

    alt(it) = radius/1000 - R_J; % km
    % solve for new alt
    if alt(it) <= alt_f % km
        fprintf('Time to decay: %.2f seconds \n',T);
        break
    end
end

figure;
subplot(2,1,1);
semilogy((1:it-1)*dT_sec,rho,'k--','LineWidth',2);
grid on;
xlabel('Time [sec]');
ylabel('\rho [kg / m^3]');
title('Atmospheric Density');
xlim([0 T]);

subplot(2,1,2);
plot((1:it)*dT_sec,alt,'k--','LineWidth',2);
grid on;
xlabel('Time [sec]');
ylabel('Altitude [km]');
title('Altitude');
xlim([0 T]);