%prob 1 v2 Earth-Moon
mu = 1.215*10^-2;

X1 = 0.813238725380;
V1 = 0.23757159412654352;

X2 = 1.108266892072;
V2 = 0.23344748256424838;

[L1,L2,L3] = getLpoints(mu);

% L1
X0 = [X1, 0, 0, 0, V1, 0];
p.mu = mu;
tspan = linspace(0,3,100);
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T1,X1] = ode45(@EOM_3body, tspan, X0, opts, p);

C1=jacobiConst_2D(X0(1:3)',X0(4:6)',mu);

% L2 
X0 = [X2, 0, 0, 0, V2, 0];
p.mu = mu;
tspan = linspace(0,4,100);
[T2,X2] = ode45(@EOM_3body, tspan, X0, opts, p);

C2=jacobiConst_2D(X0(1:3)',X0(4:6)',mu);

figure;
plot(X1(:,1),X1(:,2),'-','LineWidth',2); 
hold on; grid on;
plot(X2(:,1),X2(:,2),'-','LineWidth',2); 
plot(L1,0,'kx','LineWidth',2); % L1
plot(L2,0,'ko','LineWidth',2); % L2
title('2 Periodic Orbits of Equal Energy');
legend('Orbit about L1','Orbit about L2'); 
xlabel('X [normalized]'); ylabel('Y [normalized]');

x1 = X1(:,1); y1 = X1(:,2);
u1 = X1(:,4); v1 = X1(:,5);
x2 = X2(:,1); y2 = X2(:,2);
u2 = X2(:,4); v2 = X2(:,5);
%% Lagrange 1
% L1 interior unstable
n = length(X1); h = 2; eps = -0.01; tf = 3;
figure;  
subplot(2,2,1); hold on; grid on;
title('Interior Unstable about L1');
for i=1:h:n
   ic = [x1(i)+eps,y1(i),0,u1(i),v1(i),0]; 
   tspan = linspace(0,tf,100); 
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),X_interior(:,2),'r');
end
plot(X1(:,1),X1(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
% L1 interior stable
eps = -0.01;
subplot(2,2,2); hold on; grid on;
title('Interior Stable about L1');
for i=1:h:n
   ic = [x1(i)+eps,y1(i),0,u1(i),v1(i),0]; 
   tspan = linspace(0,tf,100);
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),-X_interior(:,2),'g');
end
plot(X1(:,1),X1(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
% L1 ext unstable
eps = 0.001;
subplot(2,2,3); hold on; grid on;
title('Exterior Unstable about L1');
for i=1:h:n
   ic = [x1(i)+eps,y1(i),0,u1(i),v1(i),0]; 
   tspan = linspace(0,tf,100); 
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),X_interior(:,2),'r');
end
plot(X1(:,1),X1(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
% L1 int stable
eps = 0.001;
subplot(2,2,4); hold on; grid on;
title('Exterior Stable about L1');
for i=1:h:n
   ic = [x1(i)+eps,y1(i),0,u1(i),v1(i),0]; 
   tspan = linspace(0,tf,100);
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),-X_interior(:,2),'g');
end
plot(X1(:,1),X1(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
%% Lagrange 2
% L2 interior unstable
eps = -0.0001; tf = 4;
figure; 
subplot(2,2,1); hold on; grid on;
title('Interior Unstable about L2');
for i=1:h:n
   ic = [x2(i)+eps,y2(i),0,u2(i),v2(i),0]; 
   tspan = linspace(0,tf,100);
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),X_interior(:,2),'r');
end
plot(X2(:,1),X2(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
% L2 interior stable
eps = -0.0001;
subplot(2,2,2);  hold on; grid on;
title('Interior Stable about L2');
for i=1:h:n
   ic = [x2(i)+eps,y2(i),0,u2(i),v2(i),0]; 
   tspan = linspace(0,tf,100);
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),-X_interior(:,2),'g');
end
plot(X2(:,1),X2(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
% L2 exterior unstable
eps = 0.0001; 
subplot(2,2,3); hold on; grid on;
title('Exterior Unstable about L2');
for i=1:h:n
   ic = [x2(i)+eps,y2(i),0,u2(i),v2(i),0]; 
   tspan = linspace(0,tf,100);
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),X_interior(:,2),'r');
end
plot(X2(:,1),X2(:,2),'k.'); axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');
% L2 exterior stable
eps = 0.0001;
subplot(2,2,4);  hold on; grid on;
title('Exterior Stable about L2');
for i=1:h:n
   ic = [x2(i)+eps,y2(i),0,u2(i),v2(i),0]; 
   tspan = linspace(0,tf,100);
   [~,X_interior] = ode45(@EOM_3body, tspan, ic, opts, p);
   plot(X_interior(:,1),-X_interior(:,2),'g');
end
plot(X2(:,1),X2(:,2),'k.'); 
axis equal
xlabel('X [normalized]'); ylabel('Y [normalized]');