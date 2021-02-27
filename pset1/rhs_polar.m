function zdot = rhs_polar(t,z,p)
r = z(1);
theta = z(2);
rdot = z(3);
thetadot = z(4);

rddot = r*thetadot^2-p.mu_AUD/r^2;
thetaddot = -2*rdot*thetadot/r;
zdot = [rdot;thetadot;rddot;thetaddot];