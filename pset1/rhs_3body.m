function zdot = rhs_3body(t,z,p)

x1 = z(1);
y1 = z(2);
x2 = z(3);
y2 = z(4);
x3 = z(5);
y3 = z(6);
den12 = ((x1-x2)^2 + (y1-y2)^2)^(3/2);
den13 = ((x1-x3)^2 + (y1-y3)^2)^(3/2);
den23 = ((x2-x3)^2 + (y2-y3)^2)^(3/2);
% body 1
xdot1 = -p.G*p.M2*(x1-x2)/den12 - p.G*p.M3*(x1-x3)/den13;
ydot1 = -p.G*p.M2*(y1-y2)/den12 - p.G*p.M3*(y1-y3)/den13;
% body 2
xdot2 = p.G*p.M1*(x1-x2)/den12 + p.G*p.M3*(x3-x2)/den23;
ydot2 = p.G*p.M1*(y1-y2)/den12 + p.G*p.M3*(y3-y2)/den23;
% body 3
xdot3 = p.G*p.M1*(x1-x3)/den13 + p.G*p.M2*(x2-x3)/den23;
ydot3 = p.G*p.M1*(y1-y3)/den13 + p.G*p.M2*(y2-y3)/den23;

zdot = [z(7:12);xdot1;ydot1;xdot2;ydot2;xdot3;ydot3];