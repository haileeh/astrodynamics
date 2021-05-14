% prob 1
% Define constants
mEarth = 5.9722*10^24; %[kg]
mSun = 1.989*10^30; %[kg]
mu = mEarth/(mSun+mEarth);

% X1 = 0.8132;
% V1 = 0.2375;
% X2 = 1.108;
% V2 = 0.233;
% 
% EE = [3.16,3.15,3.14];

%X0 = [1.011048635590035, 0, 0, 0, -0.010651851143342, 0];
X0 = [1.011220338466487, 0, 0, 0, -0.009041455692324, 0];
% Lagrange points
[L1,L2,L3] = getLpoints(mu);

% plot orbit
%ic =[X1 0 0 0 V1 0];

%mu = 0.04;
[X0,t_half]=differentialControl(mu,X0,3);

p.mu = mu;
tspan = linspace(0,2*t_half,100);
opts = odeset;
[T,X2] = ode45(@EOM_3body, tspan, X0, opts, p);

C2=jacobiConst_2D(X0(1:3)',X0(4:6)',mu)

%X0 = [0.988936913869454, 0, 0, 0, 0.009643328819185, 0];
X0 = [0.988882681939741 , 0, 0, 0, 0.008891028426809, 0];

r1=sqrt((mu+X0(1))^2);
r2=sqrt((X0(1)-(1-mu))^2);
vy = sqrt(2*(2*((X0(1)^2)/2 + (1-mu)/r1 + mu/r2) - C2));
X0(5) = vy;
[X0,t_half]=differentialControl(mu,X0,3);
tspan = linspace(0,2*t_half,100);
opts = odeset;
[T,X1] = ode45(@EOM_3body, tspan, X0, opts, p);
C1=jacobiConst_2D(X0(1:3)',X0(4:6)',mu)

figure;
plot(X2(:,1),X2(:,2)); 
hold on; grid on;
plot(X1(:,1),X1(:,2)); 
plot(L1,0,'kx','LineWidth',2); % L1
plot(L2,0,'ko','LineWidth',2); % L2

% now find stable and unstable interior/exterior manifolds 