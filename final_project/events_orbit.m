function [value,isterminal,direction] = events_orbit(t,z,p)
% Event function 
% Locate the time when the object returns closest to the initial point ic
% and starts to move away, and stop integration.  Also locate the time when
% the object is farthest from the initial point ic and starts to move closer.
%
% The current distance of the body is
%
%   DSQ = (y(1)-y0(1))^2 + (y(2)-y0(2))^2 = <y(1:2)-y0(1:2),y(1:2)-y0(1:2)>
%
% A local minimum of DSQ occurs when d/dt DSQ crosses zero heading in
% the positive direction.  We can compute d/dt DSQ as
%
%   d/dt DSQ = 2*(y(1:2)-y0)'*dy(1:2)/dt = 2*(y(1:2)-y0)'*y(3:4)


dDSQdt = 2 * ((z(1:3)-p.ic(1:3)')' * z(4:6));
value = [dDSQdt; dDSQdt];
isterminal = [1;  0];         % stop at local minimum
direction  = [1; -1];         % [local minimum, local maximum]
end