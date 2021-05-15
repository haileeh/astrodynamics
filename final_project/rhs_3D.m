function Zdot = rhs_3D(t,z,p)

r = ((z(1)^2+z(2)^2+z(3)^2))^(3/2);
xdot = -p.mu*z(1)/r;
ydot = -p.mu*z(2)/r;
zdot = -p.mu*z(3)/r;
Zdot = [z(4);z(5);z(6);xdot;ydot;zdot];