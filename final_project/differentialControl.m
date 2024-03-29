function [X0,t_half]=differentialControl(mu,X0,t_half)
%% differential correction 
% Given a guess for the initial conditions, differential correction finds
% the correct initial conditions that will produce a trajectory to the
% desired final state

% For halo orbit, I set the desired final state as the initial state
% Use symmetry about the x-z plane and periodicity

% Yields initial conditions and t_half (or just t)
% example that works
% mu = 0.04; X0 = [1.14 0 0.28 0 -0.316028 0]; t_half = 3;

p.mu = mu;%0.04; %mu; %3.986004418*10^5; %[km^3 s^-2]
it = 0;
tol = 1e-14;%1e-8;
tf = t_half;%1.433655;
while 1
    tspan = linspace(0,tf,2000);
    Phi0 = eye(6);
    ic = [X0, reshape(Phi0, 1, 6^2)];
    p.ic = ic;
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'Events', @y_cross); 
    [T,X] = ode45(@EOM_3body_var, tspan, ic, opts, p);

    x_t1 = X(end,1:6);
    del_x1 = X0 - x_t1; % desired is initial state
    
    it = it + 1;
    if it > 10 % done yet?
        disp('INCOMPLETE')
        break
    end
    
    % check tolerance
    if abs(x_t1(4)) < tol && abs(x_t1(6)) < tol
        fprintf('converged\n')
        t_half = T(end);
        break;
    else
        % continue
        Phi = reshape(X(end, 7:end), 6, 6);
        Xdot = EOM_3body_var(0,X(end,:)',p);           
         % Howell 1984 paper - equivalent to Sanchez
%          rhsMat = [Phi(4,3), Phi(4,5); Phi(6,3), Phi(6,5)] - (1/Xdot(2))*[Xdot(4); Xdot(6)]*[Phi(2,3), Phi(2,5)];
%          delta = rhsMat\[del_x1(4);del_x1(6)];
%          X0(3) = X0(3) + delta(1);
%          X0(5) = X0(5) + delta(2);

        % Sanchez
        rhsMat = [Phi(2,3), Phi(2,5) Xdot(2); Phi(4,3), Phi(4,5), Xdot(4);...
            Phi(6,3), Phi(6,5) Xdot(6)];
        if rcond(inv(rhsMat)) < 1e-16
            keyboard
        end
        deltaState = rhsMat\[0; del_x1(4); del_x1(6)];
        X0(3) = X0(3) + deltaState(1); %z0
        X0(5) = X0(5) + deltaState(2); %ydot0
        tf = tf+deltaState(3); % don't need to do this because of y_cross events
    end
end