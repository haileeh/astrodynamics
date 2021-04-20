%% prob 4 - porkchop plots
% Mars to Earth between 2021 and 2031
% Need launch dates, arrival dates, delta V associated
% plot the Hohmann transfer launch and arrival time

% constants
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
muSun_hr = muSun_s * 3600^2; %[km^3/hr^2]
muSun_day = muSun_hr * 24^2; %[km^3/day^2]
rEarth_Sun = 149.60*10^6; %[km]
rMars_Sun = 227.9*10^6; %[km]
mu_n = muSun_day / rEarth_Sun^3; %1/days^2
v_Mars = sqrt(muSun_s/rMars_Sun); %km/sec

% mean motion of Earth
r1 = rEarth_Sun/rEarth_Sun; % non-dimensional
omega1 = sqrt(mu_n/r1^3); %rad/days
T_E = 2*pi/omega1; % days
nEarth = 2*pi/T_E; %[rad/day]
nEarth_yr = 2*pi; %[rad/year]
% mean motion of Mars
r2 = rMars_Sun/rEarth_Sun; % non-dimensional
omega2 = sqrt(mu_n/r2^3); % rad/days
T_M = 2*pi/omega2; % days
nMars = 2*pi/T_M; %rad/days
nMars_yr = nMars * T_E; %rad/yrs

% find transfer time in days
aH = (rEarth_Sun+rMars_Sun)/2; % [km]
tH = 0.5*sqrt(4*pi^2*aH^3/ muSun_day);
tH_yr = tH / T_E;

%% Results for Hohmann Transfer (Mars to Earth)
% last transfer day
yr_LeavePrev = 2020.4679;
yr_ArrPrev = yr_LeavePrev + tH_yr;

% What is the synodic period? Mars orbit relative to Earth
T_syn_days = T_E*T_M / (abs(T_E-T_M)); 
T_syn_years = T_syn_days / T_E;

% Earth to Mars!
leaveDates = yr_LeavePrev + T_syn_years*(-1:1:6);
arriveDates = leaveDates + tH_yr;

% find Mars to Earth launch dates and arrival dates
% phase angle between earth and Mars when spacecraft reaches Mars is
phi_f = pi - nEarth*tH; %[rad]
N = 1;
t_wait = ( -2*phi_f - 2*pi*N ) / (nMars - nEarth); % [days]
leaveMars = arriveDates + t_wait/T_E;
arriveEarth = leaveMars + tH_yr;

%% solve for porkchop plot
v_hohmann = sqrt(muSun_s/rMars_Sun)*(1-sqrt(2*rEarth_Sun/(rEarth_Sun+rMars_Sun))); %v_mars - v_pi
maxDV = 1.5*v_hohmann;
startDates = leaveMars(1):(1/T_E):2031;
tof = 10:1:500;
b = 0; DVarray=[]; tofArray=[];
startArray=[];
for i=1:length(startDates)
    pos_M = rMars_Sun*[cos(nMars_yr*startDates(i));sin(nMars_yr*startDates(i));0];
    vel_M = v_Mars*[-sin(nMars_yr*startDates(i));cos(nMars_yr*startDates(i));0];
    for j=1:length(tof)
        stDates = startDates(i) + tof(j)/T_E;
        pos_E = rEarth_Sun*[cos(nEarth_yr*stDates);sin(nEarth_yr*stDates);0];
        [launchVel_short,vEshort] = lambert(pos_M', pos_E', tof(j), 0, muSun_s); % yields [km/s]
        [launchVel_long,vElong] = lambert(pos_M', pos_E', -tof(j), 0, muSun_s);
        DVshort = norm(launchVel_short'-vel_M);
        DVlong = norm(launchVel_long'-vel_M);
        % check against maxDV
        if min([DVshort, DVlong]) <= maxDV
            % save values for plot
            if DVshort <= maxDV
                b=b+1;
                DVarray(b) = DVshort;
                tofArray(b) = tof(j);
                startArray(b) = startDates(i);
            end
            if DVlong <= maxDV
                b=b+1;
                DVarray(b) = DVlong;
                tofArray(b) = tof(j);
                startArray(b) = startDates(i);
            end
        end
    end
end

%% find V_hohmann in data
jj = 0;
dvHo = []; 
idxDV=[];
for k=1:length(DVarray)
   if DVarray(k) <= v_hohmann*1.0001 
       if DVarray(k) >= 0.9999*v_hohmann
           jj = jj+1;
           dvHo(jj) = DVarray(k);
           idxDV(jj) = k;
       end
   end
end

%% find maxDV in data - outline
idxMax = [];
ii = 0;
for k=1:length(DVarray)
    if DVarray(k) <= 1.01*maxDV
        if DVarray(k) >= 0.99*maxDV
            ii = ii+1;
            idxMax(ii) = k;
        end
    end
end

%% plot
endArray = startArray + tofArray/T_E;
figure; 
%plot(startArray, endArray,'o');
hold on; grid on;
plot(startArray(idxMax), endArray(idxMax),'.');
plot(startArray(idxDV),endArray(idxDV),'.','LineWidth',5);
xlim([2021 2031]);
ylim([2021 2031]);
xlabel('Launch Date [yr]');
ylabel('Arrival Date [yr]');
title('Mars to Earth Porkchop Plot');
legend('1.5*\Delta V_{Hohmann}','\Delta V_{Hohmann}','Location','Southeast')