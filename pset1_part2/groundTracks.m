function groundTracks(Omega_rad,inc_rad,omega_rad,X3,Y3,t)

persistent t_Earth
if isempty(t_Earth)
    % Earth's rotational period
    t_Earth = 23*3600 + 56*60+ 4.098;
end
% rotation from perifocal to inertial space in Earth-fixed coordinate
% system
TT1 = [cos(Omega_rad), -sin(Omega_rad), 0; sin(Omega_rad), cos(Omega_rad), 0; 0 0 1];
TT2 = [1 0 0; 0 cos(inc_rad), -sin(inc_rad); 0, sin(inc_rad), cos(inc_rad)];
TT3 = [cos(omega_rad), -sin(omega_rad), 0; sin(omega_rad), cos(omega_rad), 0; 0 0 1];
T = TT1*TT2*TT3;
for i=1:length(t)
    vec = T*[X3(i);Y3(i);0]; 
    X(i) = vec(1); Y(i) = vec(2); Z(i) = vec(3);
    beta = 2*pi*t(i)/t_Earth;  % t is time history
    beta0 = 0; % specific reference location
    T0 = [cos(beta)+beta0, sin(beta)+beta0, 0;...
        -sin(beta)+beta0, cos(beta)+beta0, 0;...
        0 0 1];
    sol = T0*[X(i);Y(i);Z(i)];
    X0(i) = sol(1); Y0(i) = sol(2); Z0(i) = sol(3);
    delta(i) = (180/pi) * asin(Z0(i) / sqrt(X0(i)^2+Y0(i)^2+Z0(i)^2)); % latitude
    phi(i) = (180/pi) * atan2(Y0(i),X0(i)); % longitude
end

figure;
plot(phi, delta,'o'); hold on; grid on;
title('Ground Tracks');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
xlim([-180 180]);
ylim([-90 90]);

figure; 
subplot(1,2,1);
plot3(X,Y,Z);hold on; grid on;
plot3(X(1),Y(1),Z(1),'g*');
plot3(X(end),Y(end),Z(end),'rs');
axis square
title('Orbit fixed in space');

subplot(1,2,2);
plot3(X0,Y0,Z0);hold on; grid on;
plot3(X(1),Y(1),Z(1),'g*');
plot3(X(end),Y(end),Z(end),'rs');
axis square
title('Ref Frame Rotating Earth'); 

figure
subplot(2,1,1);
plot(t,phi); title('Longitude');
subplot(2,1,2)
plot(t,delta); title('Latitude');