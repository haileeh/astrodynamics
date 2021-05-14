function [X0_new,t_half_new] = continuation(X0_1,X0_2,mu)
% compute discrete sampling of solutions across several familiies of
% interest
% linearly extrapolate 2 converged ICs
% set fixed state (x or z) as the more rapidly changing of the 2

% constant
eps = -0.001;% TODO how to bound this
x01 = X0_1(1); x02 = X0_2(1);
z01 = X0_1(3); z02 = X0_2(3);
ydot01 = X0_1(5); ydot02 = X0_2(5);

if abs(x02-x01) > abs(z02-z01) % x is fixed state
    z0_new = (z02 - z01)/(x02-x01)*eps + z02;
    ydot0_new = (ydot02 - ydot01)/(x02-x01)*eps + ydot02;
    X0_new = [x02+eps, 0, z0_new, 0, ydot0_new, 0];
else % z is fixed state
    x0_new = (x02 - x01)/(z02-z01)*eps + x02;
    ydot0_new = (ydot02 - ydot01)/(z02-z01)*eps + ydot02;
    X0_new = [x0_new, 0, z02+eps, 0, ydot0_new, 0];
end

% call differential corrector
[X0_new,t_half_new]=differentialControl_zFixed(mu,X0_new,3); % replace initial time guess (or don't use this as input)
