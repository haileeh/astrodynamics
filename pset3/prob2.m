% Lambert's problem
rEarth_Sun = 149.60*10^6; %[km]
rMars_Sun = 227.9*10^6; %[km]
%% prob 2.1
% calculate properties of an orbit to Mars that intercepts circular orbit:
r1 = 1;
R1 = [r1; 0; 0];
r2 = 1.52; %rMars/rEarth
theta = 150; %[deg]
R2 = [r2*cosd(theta); r2*sind(theta); 0];
T1 = 3.5;
mu = 1;

zinit = 6;%6.376;
[v1,v2]=myLambert(R1,R2,theta,zinit,T1,mu);

% Step 8. orbital elements using R1 and v1
r = r1;
v = norm(v1);
v_r = dot(R1,v1)/r1;
hvec = cross(R1,v1);
h = norm(hvec);
i = acosd(hvec(3)/h);
Nvec = cross([0;0;1],hvec);
N = norm(Nvec);
Omega = acosd(Nvec(1)/N);
e_vec = 1/mu * ((v^2 - mu/r1)*R1 - r1*v_r*v1);
e = norm(e_vec);
omega = acosd(dot(Nvec,e_vec/(N*e)));
th = acosd(dot(e_vec,R1)/(e*r1));

% need to solve for semi-major axis 
C = R2-R1;
c = norm(C);
s = 0.5*(r1+r2+c);
syms smaxis real
alpha = 2*asin(sqrt(s/(2*smaxis)));
beta = 2*asin(sqrt((s-c)/(2*smaxis)));
T = smaxis^(3/2)/sqrt(mu) * (alpha - beta - (sin(alpha) - sin(beta)));
eqn = T-T1;
a = vpasolve(eqn,smaxis); %correct
%% prob 2.2 - determine T2, longer time for a transfer ellipse of same
% semi-major axis
% Determine velocities for each orbit at the inititation and at the arrival
alpha0 = 2*asin(sqrt(s/(2*a)));
alpha = 2*pi-alpha0;
beta = 2*asin(sqrt((s-c)/(2*a)));
T2f = a^(3/2)/sqrt(mu) * (alpha - beta - (sin(alpha) - sin(beta)));

[v1_f,v2_f]=myLambert(R1,R2,360-theta,zinit,T2f,mu);

%% prob 2.3 solve 2.1 using ODE
% verify correct time and correct arrival velocities
inits = [R1(1); R1(2); v1(1); v1(2)];
tspan = linspace(0,T1,100);
opts = odeset;

p.mu = 1;

[t,z] = ode45(@rhs,tspan,inits,opts,p);

x = z(:,1); y = z(:,2);
vx = z(:,3); vy = z(:,4);
% matches R2 and v2
diff_x = x(end)-R2(1)
diff_y = y(end)-R2(2)
diff_vx = vx(end)-v2(1)
diff_vy = vy(end)-v2(2)
keyboard
%% check using lambert solver - can't get this to match
R1d = R1*rEarth_Sun;
R2d = R2*rEarth_Sun;
muSun_s = 1.327*10^20 / 1000^3; %[km^3/s^2]
sqrtMuSun = sqrt(2.96*10^-4); % [1/day]
sqrtMuSun_s = sqrtMuSun * 1/(24*3600) ; 
% muSun_hr = muSun_s * 3600^2; %[km^3/hr^2]
% muSun_day = muSun_hr * 24^2; %[km^3/day^2]
% muSun_day2 = muSun_day / rEarth_Sun^3 %[1/day^2];
mu = muSun_s;%mu * rEarth_Sun^3 * (sqrtMuSun /(24*3600))^2;
mu = 1 * rEarth_Sun^3 * sqrtMuSun_s^2;
T = T1 / sqrtMuSun;
[V1, V2, extremal_distances, exitflag] = lambert(R1d', R2d', T, 0, mu);
v1-V1'
v2-V2'