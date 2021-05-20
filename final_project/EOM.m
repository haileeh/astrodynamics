% 1. Define coordinate frames and equations of motion for the two-body 
% problem in low Earth orbit and for the three-body problem as spacecraft 
% escape from the Earth’s sphere of influence towards Sun-Earth L2

%% 2-Body Problem
% need to define:Omega_rad,inc_rad,omega_rad,X3,Y3,t

% Perifocal to ECEF to ECI
% Earth's rotational period
t_Earth = 23*3600 + 56*60+ 4.098;
TT1 = [cos(Omega_rad), -sin(Omega_rad), 0; sin(Omega_rad), cos(Omega_rad), 0; 0 0 1];
TT2 = [1 0 0; 0 cos(inc_rad), -sin(inc_rad); 0, sin(inc_rad), cos(inc_rad)];
TT3 = [cos(omega_rad), -sin(omega_rad), 0; sin(omega_rad), cos(omega_rad), 0; 0 0 1];
T = TT1*TT2*TT3;

for i=1:length(t)
    vec = T*[X3(i);Y3(i);0]; 
    X(i) = vec(1); Y(i) = vec(2); Z(i) = vec(3); % ECEF
    beta = 2*pi*t(i)/t_Earth;  % t is time history
    beta0 = 0; % specific reference location
    T0 = [cos(beta)+beta0, sin(beta)+beta0, 0;...
        -sin(beta)+beta0, cos(beta)+beta0, 0;...
        0 0 1];
    sol = T0*[X(i);Y(i);Z(i)];
    X0(i) = sol(1); Y0(i) = sol(2); Z0(i) = sol(3); % ECI
end



%% 3-Body Problem
% Heliocentric, sun at origin O, inertial
% ihat points towards First Point of Aries during specific epoch on Jan
% 2000
% have rhs3D function (1-body), can make 2-body and 3-body equivalents