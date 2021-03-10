%% prob 1.1 - orbit determination
% constants from problem
R0 = [0.0468388; 1.60929; 0];
V0 = [-0.0134275; 0.00898728; 0]; %first dim should be positive
mu = 0.00295;
V_norm = norm(V0);
R_norm =  norm(R0);

% eccentricity vector
e_vec = cross(V0,cross(R0,V0))/mu - R0/R_norm;
e = norm(e_vec);

% find semi-major axis
a = (2 / R_norm - V_norm^2/mu)^-1;

% find angular momentum
h_vec = cross(R0,V0);
h = norm(h_vec);

% parameter of ellipse
p_param = a*(1-e^2);

%% prob 1.2 - coordinate system
ihatE = e_vec/e;
ihatH = h_vec/h;
ihatP = cross(ihatE,ihatH);

A = [ihatE, ihatP, ihatH]; 
r0 = A*R0;
v0 = A*V0;

%% prob 1.3 - what is E0?
sigma0 = dot(R0,V0) / sqrt(mu);
r0_norm = norm(r0);
v0_norm = norm(v0);
sigma = dot(r0,v0) / sqrt(mu);
E0 = atan2((sigma0/sqrt(a)),(1-R_norm/a));%atan((sigma0/sqrt(a))/(1-R_norm/a));%atan2(1-R_norm/a, sigma0/sqrt(a));

%% prob 1.4 - express values of x, y coordinate using E
expf02 = sqrt((1+e) / (1-e)) * tan(E0/2);
f0 = 2 * atan(expf02);

n = sqrt(mu/a^3);
t0 = (1/n) * (E0 - e*sin(E0));
%syms Et real
Et = linspace(E0,E0+9,100)';
Etilde = Et-E0;
%func = (1/n)*(t-t0) - Etilde.*(1-r0_norm/a).*sin(Etilde) + (sigma0/sqrt(a)).*(cos(Etilde)-1);
r = a + (r0_norm-a)*cos(Etilde) + sigma0*sqrt(a)*sin(Etilde);
t = (Et - e*sin(Et))/n;
F = 1 - (a/r0_norm)*(1-cos(Etilde));
G = (t-t0) + sqrt(a^3/mu)*(sin(Etilde)-Etilde);
Fdot = -sqrt(mu*a)/(r*r0_norm)*sin(Etilde);
Gdot = 1 - (a/r)*(1-cos(Etilde));
x0 = r0(1);
y0 = r0(2);
xdot0 = v0(1);
ydot0 = v0(2);
x = x0*F + xdot0*G;
y = y0*F + ydot0*G;

figure; 
plot(x,y); hold on;
xlabel('X-Position');
ylabel('Y-Position');
xlim([-1.8,0.2]); ylim([-0.9,1.1]);
axh = gca; % use current axes
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', 'k', 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', 'k', 'LineStyle', linestyle);

%% prob 1.5 - perform an initial value calculation
% plot the orbit
p.mu_AUD = mu;
inits = [x0 y0 xdot0 ydot0];
% time span (one period)
tf = 2*pi*a^(3/2)/sqrt(p.mu_AUD);
tspan = linspace(0,tf,300);
% solve ode
opts = odeset;
[~,z]=ode45(@rhs,tspan,inits,opts,p);

figure
plot(z(:,1),z(:,2)); hold on;
title('X-Position vs. Y-Position');
xlim([-1.8,0.2]); ylim([-0.9,1.1]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', 'k', 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', 'k', 'LineStyle', linestyle);

% plot on top of each other?

%% prob 1.6 - series expansion
% doesn't quite look right
clear t
%tf = 132;
t = linspace(t0,4*tf,100); %should I use a different tf to make orbit complete?
h0 = mu/r0_norm^3;
p0 = dot(r0,v0)/r0_norm^2;
f = 1 - h0*((t - t0).^2)/2 + h0*p0*((t - t0).^3)/2;
g = (t - t0) - h0*((t - t0).^3)/6;
x2 = f*x0 + g*xdot0;
y2 = f*y0 + g*ydot0;

figure
plot(x2,y2); hold on;
title('X-Position vs. Y-Position');
xlim([-1.8,0.2]); ylim([-0.9,1.1]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', 'k', 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', 'k', 'LineStyle', linestyle);