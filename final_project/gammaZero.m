function [value,isterminal,direction] = gammaZero(t,z,p)

x1 = z(1); y1=z(2); z1=z(3);
xdot1 = z(4); ydot1 = z(5); zdot1 = z(6);
mu = p.mu;
R=p.R;
ha = sqrt( (x1-1+mu)^2 + y1^2 + z1^2) - R;
gamma = acos( sqrt( (y1*zdot1 + z1*ydot1)^2 + ...
    (z1*xdot1 - (x1-1+mu)*zdot1)^2*( (x1-1+mu)*ydot1 - y1*xdot1)^2) /...
        ( (ha + R)*sqrt(xdot1^2+ydot1^2+zdot1^2)) );
value = gamma;
isterminal = 1;         % stop at local minimum
direction  = 1;         % [local minimum, local maximum]