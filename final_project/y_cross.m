function [value,isterminal,direction] = y_cross(t,z,p)


value = [z(2);z(2)];
isterminal = [1;  1];         % stop at local minimum
direction  = [1; -1];         % [local minimum, local maximum]