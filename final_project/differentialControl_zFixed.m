function [X0,t_half]=differentialControl_zFixed(mu,X0,t_half)
%% differential correction 
% Given a guess for the initial conditions, differential correction finds
% the correct initial conditions that will produce a trajectory to the
% desired final state

% For halo orbit, I set the desired final state as the initial state
% Use symmetry about the x-z plane and periodicity

% Yields initial conditions and t_half (or just t)
% example that works

p.mu = mu;
it = 0;
tol = 1e-8;
tf = t_half;
while 1
    tspan = linspace(0,tf,2000);
    Phi0 = eye(6);
    ic = [X0, reshape(Phi0, 1, 6^2)];
    p.ic = ic;
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events', @y_cross); 
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

        % Sanchez
        rhsMat = [Phi(2,1), Phi(2,5) Xdot(2); Phi(4,1), Phi(4,5), Xdot(4);...
            Phi(6,1), Phi(6,5) Xdot(6)];
        if rcond(inv(rhsMat)) < 1e-16
            keyboard
        end
        deltaState = rhsMat\[0; del_x1(4); del_x1(6)];
        X0(1) = X0(1) + deltaState(1); %x0
        X0(5) = X0(5) + deltaState(2); %ydot0
        tf = tf+deltaState(3);
    end
end