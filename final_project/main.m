plotFlag = 0;
unstableFlag = 0;
% Define constants
mEarth = 5.9722*10^24; %[kg]
mSun = 1.989*10^30; %[kg]
mu = mEarth/(mSun+mEarth);
au = 149.5978714*10^6; %km - distance between Sun and Earth

G = (6.67408e-11); % m^3 / kg / s^2
n = -sqrt((mSun+mEarth)*G/((au)*10^3)^3); % rad/s

% amplitude of halo orbit in km
Az(1) = 125000; 
Az(2) = 150000;

% Lagrange points
[L1,L2,L3] = getLpoints(mu);

% which Lagrange pt?
L = L2;

X0 = zeros(length(Az),6);
for i=1:length(Az)
    % Third-order approximation of initial conditions
    x_richardson = richardson(mu,Az(i),n,au,L);
    
    % Converting from dimensioned coordinates about the Lagrange point
    % to be about the center of mass of the system and ND
    x_3a = [L*au+x_richardson(1); 0; x_richardson(3); 0; L*au-(L*au*n+x_richardson(5))/n;0]/au;
    
    % fine-tune initial conditions of halo orbit using differential correction
    [X0(i,:),t_half(i)]=differentialControl_zFixed(mu,x_3a',3); % TODO initial time
end
%% continuation!
for k=1:5
    [X0(2+k,:),t_half(2+k)] = continuation(X0(k,:),X0(k+1,:),mu);
end
%% plot halo orbits
if plotFlag
    figure;
    plot3(1-mu,0,0,'bo','LineWidth',2); % Earth
    hold on; grid on;
    plot3(L1,0,0,'mx','LineWidth',2); % L1
    plot3(L2,0,0,'ms','LineWidth',2); % L2
    for j=1:size(X0,1)
        ic = X0(j,1:6);
        p.mu = mu;
        tspan = linspace(0,2*t_half(j),100);
        opts2 = odeset;
        [T,X] = ode45(@EOM_3body, tspan, ic, opts2, p);
        
        plot3(X(:,1),X(:,2),X(:,3),'k');
    end
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Halo orbit family about L2');
    legend('Earth','L1','L2','Halo orbits');
    % rotates about lagrange point from y-z projection
end
%% find monodromy matrix
for j=1:1%12
    ic = [X0(j,1:6), reshape(eye(6),6^2,1)'];
    tspan = linspace(0,2*t_half(1),100);
    [T,X] = ode45(@EOM_3body_var, tspan, ic, opts2, p);
    for i=1:length(T)
        Phi(:,:,i) = reshape(X(i,7:end), 6, 6);
    end
    Phif = Phi(:,:,end);
    
    [V,D]=eig(Phif); % D is eigenvalues, V is eigenvectors (columns)
    eigD = diag(D);
    Y_u = V(:,1); %unstable
    Y_s = V(:,2); %stable
    
    %% find stable manifold
    ic_t = X(:,1:6);
    [state_tS,C_man,tMan] = findManifold(ic_t,p,Y_s',eigD(2));
    %C_man same energy for invariant manifold
    if plotFlag
        figure;
        plot3(X(:,1),X(:,2),X(:,3),'*');hold on; grid on;
        plot3(1-mu,0,0,'bo','LineWidth',2); % Earth
        plot3(L1,0,0,'mx','LineWidth',2); % L1
        plot3(L2,0,0,'kx','LineWidth',2); % L2
        for i=1:2:100
            plot3(state_tS(:,1,i),state_tS(:,2,i),state_tS(:,3,i),'g');
        end
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title('Stable Manifold');
        legend('Halo orbit','Earth','L1','L2','Stable Manifold');
    end
    %% find unstable manifold
    if unstableFlag
        ic_t = X(:,1:6);
        [state_tU,C_man,tMan] = findManifold(ic_t,p,Y_u',eigD(1));
        %C_man same energy for invariant manifold
        if plotFlag
            figure;
            plot3(X(:,1),X(:,2),X(:,3),'*'); hold on; grid on;
            plot3(1-mu,0,0,'bo','LineWidth',2); % Earth
            plot3(L1,0,0,'mx','LineWidth',2); % L1
            plot3(L2,0,0,'kx','LineWidth',2); % L2
            
            for i=1:2:100
                plot3(state_tU(:,1,i),state_tU(:,2,i),state_tU(:,3,i),'r'); 
            end
            
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Unstable Manifold');
            legend('Halo orbit','Earth','L1','L2','Unstable Manifold');
        end
    end
    %% closest to Earth
    diffSmall = 1e8;
    idx = [0 0];
    for i=1:size(state_tS,3)
        xM = state_tS(:,1,i); yM = state_tS(:,2,i); zM = state_tS(:,3,i);
        for k=1:length(xM)
            diff = norm([1-xM(k);yM(k);zM(k)]);
            if diff < diffSmall
                diffSmall = diff;
                idx = [i,k];
            end
        end
    end
    minDist = diffSmall*au;
    % closest manifold point to bring sc to halo orbit in normalized cr3bp
    % coordinates
    % state_tS dimensions = [time step, state pos, ic #]
    tti_state = state_tS(idx(2),:,idx(1)); % TTI ; X1 - is this correct?
    hoi_state = ic_t(idx(1),:); % HOI - corresponds to pt above; X0 (correct)
    hoi0_state = state_tS(1,:,idx(1));
    tf = tMan(idx(2),idx(1));
    % the velocities in hoi_state would be subscript h
    % the velocities in the corresponding pt on the manifold would be 0
    HOI_diffCorrector(tti_state,hoi_state,hoi0_state,mu,tf)
end