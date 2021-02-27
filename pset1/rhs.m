function zdot = rhs(t,z,p)

xdot = -p.mu_AUD*z(1)/((z(1)^2+z(2)^2))^(3/2);
ydot = -p.mu_AUD*z(2)/((z(1)^2+z(2)^2))^(3/2);
zdot = [z(3);z(4);xdot;ydot];