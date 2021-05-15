% Define constants
mEarth = 5.9722*10^24; %[kg]
mSun = 1.989*10^30; %[kg]
mu = mEarth/(mSun+mEarth);
% mu = 3.040357143*10^-6; % from Richardson; mu = G*E in this context; not sure what E is
au = 149.5978714*10^6; %km %R12 - distance between Sun and Earth


%n = 1.990986606*10^-7; % rad/s
G = (6.67408e-11); % m^3 / kg / s^2
n = -sqrt((mSun+mEarth)*G/((au)*10^3)^3); % rad/s

% amplitude of halo orbit
%Az(1) = 270000; %110000; %km % from Richardson's paper? nominal orbit amplitude
Az(1) = 125000; %test
Az(2) = 100000;
%Az(2) = Az(1)*1.1;

% Lagrange points
[L1,L2,L3] = getLpoints(mu);

% which Lagrange pt?
L = L2;

X0 = zeros(length(Az),6);
for i=1:length(Az)
    % Third-order approximation of initial conditions
    x_richardson = richardson(mu,Az(i),n,au,L);
    
    % convert: Converting from dimensioned coordinates about the Lagrange point
    % to be about the center of mass of the system and ND
    x_3a = [L*au+x_richardson(1); 0; x_richardson(3); 0; L*au-(L*au*n+x_richardson(5))/n;0]/au;
    
    % fine-tune initial conditions of halo orbit using differential correction
    [X0(i,:),t_half(i)]=differentialControl_zFixed(mu,x_3a',3); % replace initial time guess (or don't use this as input)
end
%% continuation!
% loop over this
% for k=1:10
%     [X0(2+k,:),t_half(2+k)] = continuation(X0(k,:),X0(k+1,:),mu);
% end
%% plot halo orbits
figure;
for j=1:size(X0,1)
    ic = X0(j,1:6);
    p.mu = mu;
    tspan = linspace(0,2*t_half(j),100);
    opts2 = odeset;
    [T,X] = ode45(@EOM_3body, tspan, ic, opts2, p);
    
    plot3(X(:,1),X(:,2),X(:,3));
    hold on; grid on;
    %plot3(X0(j,1),X0(j,2),X0(j,3),'rx','LineWidth',2);
end
plot3(1-mu,0,0,'bo','LineWidth',2); % Earth
%plot3(-mu,0,0,'yo','LineWidth',2); % Sun
plot3(L1,0,0,'kx','LineWidth',2); % L1
plot3(L2,0,0,'ko','LineWidth',2); % L2
xlabel('X'); ylabel('Y'); zlabel('Z');
% looks like it rotates about lagrange point from y-z projection
%% find monodromy matrix of last orbit
ic = [ic, reshape(eye(6),6^2,1)'];
[T,X] = ode45(@EOM_3body_var, tspan, ic, opts2, p);
for i=1:length(T)
    Phi(:,:,i) = reshape(X(i,7:end), 6, 6);
end
Phif = Phi(:,:,end);
[V,D]=eig(Phif); % D is eigenvalues, V is eigenvectors (columns)
Y_u = V(:,1); %unstable
Y_s = V(:,2); %stable

%% find stable manifold
extFlag = 0; % remove stable flag, replace with eig vector
ic_t = X(:,1:6); % would need to save off more trajectories; or re-solve from one of the ICs (X0)
[state_t,C_man] = findManifold(ic_t,extFlag,p,Y_s');

figure; 
plot3(X(:,1),X(:,2),X(:,3),'*');
hold on; grid on;
for i=1:2:100 %50 to 100 for only inward
   plot3(state_t(1,1,i),state_t(1,2,i),state_t(1,3,i),'ro');
   plot3(state_t(:,1,i),state_t(:,2,i),state_t(:,3,i),'g');  % green for stable
end
plot3(L1,0,0,'kx','LineWidth',2); % L1
plot3(L2,0,0,'ko','LineWidth',2); % L2
xlabel('X'); ylabel('Y'); zlabel('Z');
plot3(1-mu,0,0,'bo','LineWidth',2); % Earth
%plot3(-mu,0,0,'yo','LineWidth',2); % Sun

%% interception with parking orbit
rEarth = 6378.1; %km
pOrbit.mu = 3.986004418*10^5/rEarth^3;% / au^3; %1/s^2, non-dim wrt au
pOrbit.r0 = [35000,0,0]/rEarth;%/au; %1000 is good
v = sqrt(pOrbit.mu/norm(pOrbit.r0));
pOrbit.v0 = [0,v,0];%/au; %1/sec
[parkingOrbX,parkingOrbY,parkingOrbZ] = closestPoint(pOrbit,state_t,au,rEarth);
plot3(parkingOrbX,parkingOrbY,parkingOrbZ,'yo');