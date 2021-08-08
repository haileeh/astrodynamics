% from: https://gist.github.com/ChristopherRabotin/8293e36daf5454774ce34b5103842273
function state = halofunc(X)
% Parse state vector into variables
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);
Phi = reshape(X(7:end), 6,6);

% Simplication
mu = 3.986004418*10^5; %[km^3 s^-2]
r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x+mu-1)^2 + y^2 + z^2);

ax = x + 2*vy - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
ay = y - 2*vx - (1-mu)*y/r1^3 - mu*y/r2^3;
az = -(1-mu)*z/r1^3 - mu*z/r2^3;
% Compute STM
Delta = zeros(3,3); Delta(1,2) = 1; Delta(2,1) = -1;
fancyU = [Omegaxx, Omegaxy, Omegaxz; Omegayx, Omegayy, Omegayz; Omegazx, Omegazy, Omegazz];
A_matrix = [zeros(3,3), eye(3); -fancyU, 2*Delta];
A = A_matrix(X(1:6)', mu);
Phid = A*Phi;

state = [vx; vy; vz; ax; ay; az; reshape(Phid, 6^2, 1)];

end

function A = A_matrix(X6,mu)
A = zeros(6,6);
Delta = zeros(3,3); Delta(1,2) = 1; Delta(2,1) = -1;

x = X6(1);
y = X6(2);
z = X6(3);
x2 = x^2;
y2 = y^2;
z2 = z^2;
r2 = x2 + y2 + z2;
r232 = r2^(3/2);
r252 = r2^(5/2);

% Add the body perturbations - is this accurate in this case
dAxDx = 3 * mu * x2/r252 - mu / r232;
dAxDy = 3 * mu * x * y / r252;
dAxDz = 3 * mu * x * z / r252;
dAyDx = 3 * mu * x * y / r252;
dAyDy = 3*mu*y2/r252 - mu/r232;
dAyDz = 3 * mu * y * z / r252;
dAzDx = 3 * mu * x * z / r252;
dAzDy = 3 * mu * y * z / r252;
dAzDz = 3*mu*z2/r252 - mu/r232;

A_matrix = [zeros(3,3), eye(3); -fancyU, 2*Delta];
end