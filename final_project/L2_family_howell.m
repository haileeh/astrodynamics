% create L2 family at mu = 0.04
% https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Journals/1984_CM_How.pdf
% From "Three-Dimensional, Periodic 'Halo' Orbits,
% Kathleen Connor Howell, Celestial Mechanics 32 (1984) 53-71 

% starting points - symmetric about x-z plane
x0     =   [1.057222, 1.092791, 1.140216, 1.173414, 1.220839, 1.258203];
y0     = zeros(1,6);
z0     =   [0.300720, 0.309254, 0.298898, 0.272900, 0.200987, 0.05000];
xdot0  = zeros(1,6);
ydot0  =   [-0.238026, -0.281140, -0.316028, -0.324710, -0.310434, -0.250410];
zdot0  = zeros(1,6);
Thalf0s = [1.019032, 1.205930, 1.433655, 1.562199, 1.700458, 1.791154];
T0s     = 2*Thalf0s;
nu1s    = [-0.01038,0.61156, 11.54674,35.36097,143.9507,458.2081];
nu2s    = [-1.43755,-0.71170,-0.98759,-0.61975,0.38028,0.98301];

p.mu = 0.04;
[L1,L2,L3] = getLpoints(p.mu);

opts = odeset('RelTol',1e-5,'AbsTol',1e-4);%,'Events',@y_cross);
figure;
for j=1:6
    ic = [x0(j), y0(j), z0(j), xdot0(j), ydot0(j), zdot0(j)];
    p.ic = ic;
    tspan = [0 2*Thalf0s(j)];
    [T,X] = ode45(@EOM_3body, tspan, ic, opts, p);
    plot3(X(:,1),X(:,2),X(:,3)); hold on;
end
grid on;
plot3(L2,0,0,'ko');
%plot3(-p.mu,0,0,'go');
%plot3(1-p.mu,0,0,'bo');
xlabel('X'); ylabel('Y'); zlabel('Z');