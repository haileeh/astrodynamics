%% prob 2 a
% Cartesian coordinates
% initial conditions
R_AU = 1;
V_AU = 1; 

x0 = R_AU;
y0 = 0;
xdot0 = 0;
ydot0 = 1.3*V_AU;

p.mu_AUD = 1;

a = p.mu_AUD*R_AU/(2*p.mu_AUD-R_AU*(1.3*V_AU)^2);
inits = [x0 y0 xdot0 ydot0];
tf = 2*pi*a^(3/2)/sqrt(p.mu_AUD);
tspan = linspace(0,2*tf,300);

opts = odeset;

[t,z]=ode45(@rhs,tspan,inits,opts,p);

figure
plot(z(:,1),z(:,2));
hold on;
title('X-Position vs. Y-Position from Cartesian');
xlim([-6, 2]);
ylim([-3, 3]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);

% polar coordinates
r0 = R_AU;
theta0 = 0;
rdot0 = 0;
thetadot0 = 1.3*V_AU/R_AU;

inits = [r0 theta0 rdot0 thetadot0];

[t,z]=ode45(@rhs_polar,tspan,inits,opts,p);

r = z(:,1);
theta = z(:,2);
x = zeros(length(tspan),1);
y = x;
for i=1:length(tspan)
    x(i) = r(i)*cos(theta(i));
    y(i) = r(i)*sin(theta(i));
    if theta(i) > 2*pi
        theta(i) = theta(i) - 2*pi;
        if theta(i) > 2*pi
            theta(i) = theta(i) - 2*pi;
        end
    end
end

figure
plot(x,y);
hold on;
title('X-Position vs. Y-Position from Polar');
xlim([-6, 2]);
ylim([-3, 3]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);

figure;
subplot(2,1,1)
plot(t,r);
title('Time vs. r');
xlabel('t');
ylabel('r(t)');
grid on;

subplot(2,1,2)
plot(t,theta);
title('Time vs. \theta');
xlabel('t');
ylabel('\theta(t)');
grid on;

%% prob 2 b
mu = 0.5;
r0 = [-1;-2;0];
v0 = [-0.5;0.1;0];

v_norm = norm(v0);
r_norm =  norm(r0);

e_vec = cross(v0,cross(r0,v0))/mu - r0/r_norm;

e = norm(e_vec);

% plot the orbit
p.mu_AUD = mu;
inits = [r0(1) r0(2) v0(1) v0(2)];

% define semi-major axis
a = (2 / r_norm - v_norm^2/mu)^-1;
% slope
m = e_vec(2)/e_vec(1);
xrange = -2.59:.01:1.38;
yrange = m*xrange;
% find center point for semi-major axis
syms xk real
sympref('FloatingPointOutput',true);
eqn = a - sqrt((yrange(1)-m*xk)^2 + (xrange(1)-xk)^2);
sol = solve(eqn,xk);
xk_sol = sol(2); % in range
yk_sol = m*xk_sol;

% time span
tf = 2*pi*a^(3/2)/sqrt(p.mu_AUD);
tspan = linspace(0,tf,300);
% solve ode
[t,z]=ode45(@rhs,tspan,inits,opts,p);

figure
plot(z(:,1),z(:,2));
hold on;
plot(xrange,yrange,'--k')
plot([xrange(1) xk_sol xrange(end)],[yrange(1) yk_sol yrange(end)],'x','LineWidth',3)
quiver(0,0, e_vec(1), e_vec(2),5,'LineWidth',3);
title('X-Position vs. Y-Position');
xlim([-4, 3]);
ylim([-2.5, 3.5]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);

%% prob 2 c - hyperbola
mu = 0.5;
r0 = [-1;-2;0];
v0 = [0.65;0.65;0];

v_norm = norm(v0);
r_norm =  norm(r0);

e_vec = cross(v0,cross(r0,v0))/mu - r0/r_norm;
e = norm(e_vec);

% plot the orbit
p.mu_AUD = mu;
inits = [r0(1) r0(2) v0(1) v0(2)];


% time span
tf = 100;%2*pi*a^(3/2)/sqrt(p.mu_AUD);
tspan = linspace(0,tf,300);
% solve ode
[t,z]=ode45(@rhs,tspan,inits,opts,p);

% plot velocities to find v_inf
figure;
plot(t,z(:,3)); hold on;
plot(t,z(:,4));
title('Time histories of velocities');
legend('x-velocity','y-velocity')
xlabel('Time');
ylabel('Velocity');

v_inf = sqrt(z(end,3)^2*z(end,4)^2);
% define semi-major axis
a = -mu/v_inf^2 

figure
plot(z(:,1),z(:,2));
hold on;
%plot(xrange,yrange,'--k')
%plot([xrange(1) xk_sol xrange(end)],[yrange(1) yk_sol yrange(end)],'x','LineWidth',3)
%quiver(0,0, e_vec(1), e_vec(2),5,'LineWidth',3);
title('X-Position vs. Y-Position');
%xlim([-4, 3]);
%ylim([-2.5, 3.5]);
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
