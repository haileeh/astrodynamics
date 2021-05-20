function [state_t,C] = findManifold(ic_t,p,Y_s,eigD)
% ic_t of size n x 6
n = length(ic_t);

if Y_s(1) > 0
    % neg eps
    eps = -0.0001;
else
    % pos eps
    eps = 0.0001;
end

tf = 2.5;
%tf = 10;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
if eigD < 1
    % use reverse direction
    tspan = linspace(tf,0,1000);
else
    tspan = linspace(0,tf,1000);
end
state_t = zeros(1000,6,100);
for i=1:n
   ic = ic_t(i,:) + eps*Y_s;
   [~,state_t(:,:,i)] = ode45(@EOM_3body, tspan, ic, opts, p);
end
C = zeros(n,1);
for k=1:n
    C(k)=jacobiConst(state_t(1,1:3,k),state_t(1,4:6,k),p.mu);
end