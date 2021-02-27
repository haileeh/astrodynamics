function zdot = rhs_2body(t,z,p)

den = ((z(1)-z(3))^2 + (z(2)-z(4))^2)^(3/2);
% earth
xdotE = -p.mu_Earth*(z(1)-z(3))/den;
ydotE = -p.mu_Earth*(z(2)-z(4))/den;
% moon
xdotM = p.mu_Moon*(z(1)-z(3))/den;
ydotM = p.mu_Moon*(z(2)-z(4))/den;

zdot = [z(5);z(6);z(7);z(8);xdotE;ydotE;xdotM;ydotM];