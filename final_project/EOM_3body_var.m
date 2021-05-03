function state = EOM_3body_var(t,X,p)
x = X(1);
y = X(2);
z = X(3);

mu = p.mu;

% from EOM_3body
Xdot = EOM_3body(t,X(1:6),p);

Phi = reshape(X(7:end), 6,6);
% state transition matrix
A = stm_jacobian(x,y,z,mu);
Phid = A*Phi;

state = [Xdot; reshape(Phid, 6^2, 1)];