function [state_t,C] = findManifold(ic_t,extFlag,p,Y_s)
% ic_t of size n x 6
n = length(ic_t);
x = ic_t(:,1); y = ic_t(:,2); z = ic_t(:,3);
vx = ic_t(:,1); vy = ic_t(:,2); vz = ic_t(:,3);

if extFlag
    % neg eps
    eps = -0.001;
else %interior
    % pos eps
    eps = 0.001;
end

%tf = 0.1;
tf = 2.5;
%tf = 10;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
state_t = zeros(1000,6,100);
for i=1:n
   %ic = [x(i)+eps,y(i),z(i),vx(i),vy(i),vz(i)]; 
   ic = ic_t(i,:) + eps*Y_s;
   tspan = linspace(0,tf,1000); 
   [~,state_t(:,:,i)] = ode45(@EOM_3body, tspan, ic, opts, p);
end
C = zeros(n,1);
for k=1:n
    C(k)=jacobiConst(state_t(1,1:3,k)',state_t(1,4:6,k)',p.mu);
end
% if stableFlag
%     % find stable manifold - neg y
%     state_t(:,2,:) = -state_t(:,2,:);
% else
%     % find unstable manifold - pos y
%     % do nothing
% end