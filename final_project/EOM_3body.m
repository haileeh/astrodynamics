function Xdot = EOM_3body(t,X,p)
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);

mu = p.mu;

r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x+mu-1)^2 + y^2 + z^2);
Omega_x = x - (mu*(mu + x - 1))/r2^3 + ((mu + x)*(mu - 1))/r1^3;
Omega_y = y + (y*(mu - 1))/r1^3 - (mu*y)/r2^3;
Omega_z = (z*(mu - 1))/r1^3 - (mu*z)/r2^3;

ax = 2*vy + Omega_x;
ay = -2*vx + Omega_y;
az = Omega_z;

Xdot = [vx; vy; vz; ax; ay; az];