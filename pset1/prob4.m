%% prob 4.1

M1 = 4;
x1init = 0; 
y1init = 0;
M2 = 5;
x2init = 3;
y2init = 0;
M3 = 3;
x3init = 3;
y3init = 4;
G = 1;
T = 75;

p.G = G;
p.M1 = M1;
p.M2 = M2;
p.M3 = M3;

inits = [x1init,y1init,x2init,y2init,x3init,y3init,zeros(1,6)];

tspan = linspace(0,T,1000); 

% solve ode
opts = odeset('RelTol',2e-13);
[t,z]=ode45(@rhs_3body,tspan,inits,opts,p);

x1 = z(:,1);
y1 = z(:,2);
x2 = z(:,3);
y2 = z(:,4);
x3 = z(:,5);
y3 = z(:,6);
vx1 = z(:,7);
vy1 = z(:,8);
vx2 = z(:,9);
vy2 = z(:,10);
vx3 = z(:,11);
vy3 = z(:,12);

figure;
plot(x1,y1,'LineWidth',2);
hold on;
plot(x2,y2,'LineWidth',2);
plot(x3,y3,'LineWidth',2);
title('X-Position vs. Y-Position');
xlabel('X-Position');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
legend('M1','M2','M3')

%% prob 4.2 - COM, stays in same position

MTotal = M1 + M2 + M3;
xCMInit = (M1*x1init + M2*x2init + M3*x3init)/MTotal;
yCMInit = (M1*y1init + M2*y2init + M3*y3init)/MTotal;

inits = [xCMInit,yCMInit,0,0];

xCM = zeros(length(tspan),1);
yCM = xCM;
vxCM = xCM;
vyCM = yCM;
for i=1:length(tspan)
    xCM(i) = (M1*x1(i) + M2*x2(i) + M3*x3(i))/MTotal;
    yCM(i) = (M1*y1(i) + M2*y2(i) + M3*y3(i))/MTotal;
    vxCM(i) = (M1*vx1(i) + M2*vx2(i) + M3*vx3(i))/MTotal;
    vyCM(i) = (M1*vy1(i) + M2*vy2(i) + M3*vy3(i))/MTotal;
end

% p.mu_AUD = G*MTotal; %?
% 
% tspan = linspace(0,T,1000); 
% 
% % solve ode
% opts = odeset('RelTol',2e-16);
% [t,zCM]=ode45(@rhs,tspan,inits,opts,p);
% 
% xCM = zCM(:,1);
% yCM = zCM(:,2);
% 
% vxCM = zCM(:,3);
% vyCM = zCM(:,4);

% plot position
figure;
subplot(2,1,1)
plot(tspan,xCM,'LineWidth',2);
title('COM X-Position Time History');
xlabel('Time');
ylabel('X-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
xlim([0 T])
subplot(2,1,2)
plot(tspan,yCM,'LineWidth',2);
title('COM Y-Position Time History');
xlabel('Time');
ylabel('Y-Position');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
xlim([0 T])

% plot velocity
figure;
subplot(2,1,1)
plot(tspan,vxCM,'LineWidth',2);
title('COM X-Velocity Time History');
xlabel('Time');
ylabel('X-Velocity');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
xlim([0 T])
subplot(2,1,2)
plot(tspan,vyCM,'LineWidth',2);
title('COM Y-Velocity Time History');
xlabel('Time');
ylabel('Y-Velocity');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
xlim([0 T])

%% prob 4.3
h1 = zeros(length(tspan),3);
h2 = h1;
h3 = h1;
for i=1:length(tspan)
    h1(i,:) = cross([x1(i);y1(i);0],[vx1(i);vy1(i);0]);
    h2(i,:) = cross([x2(i);y2(i);0],[vx2(i);vy2(i);0]);
    h3(i,:) = cross([x3(i);y3(i);0],[vx3(i);vy3(i);0]);
    hCM(i,:) = cross([xCM(i);yCM(i);0],[vxCM(i);vyCM(i);0]);
end
h1 = h1(:,3);
h2 = h2(:,3);
h3 = h3(:,3);

figure;
plot(tspan,h1);
hold on;
plot(tspan,h2);
plot(tspan,h3);
title('Magnitude of Angular Momentum');
xlabel('Time');
ylabel('Ang Momentum');
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
xlim([0 T])
legend('M1','M2','M3')