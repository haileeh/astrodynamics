function [x,y,z]=closestPoint(pOrbit,manifolds,au,rEarth)
% find closest point on stable manifold to Earth parking orbit

% first, integrate desired parking orbit
ic = [pOrbit.r0, pOrbit.v0];
p.mu = pOrbit.mu;
tspan = linspace(0,100,1000); % change this
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

[T,X]=ode45(@rhs_3D,tspan,ic,opts,p); %origin is Earth
% X is normalized to rEarth
X = X*rEarth;
x = X(:,1)/au + 1; y = X(:,2)/au; z = X(:,3)/au; %add 1 au
%figure; plot3(x,y,z)
% find shortest distance compared to manifold time histories
diffSmall = 1e6;
for i=1:size(manifolds,3)
    xM = manifolds(:,1,i); yM = manifolds(:,2,i); zM = manifolds(:,3,i);
    for j=1:length(x)
        for k=1:length(xM)
            diff = norm([x(j)-xM(k);y(j)-yM(k);z(j)-zM(k)]);
            if diff < diffSmall
                diffSmall = diff;
            end
        end
    end
end
diffSmall*au