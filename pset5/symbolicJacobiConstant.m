% symbolic jacobi constant
%sympref('FloatingPointOutput',false)

syms x vy real
mEarth = 5.9722*10^24; %[kg]
mSun = 1.989*10^30; %[kg]
mu = mEarth/(mSun+mEarth);
%the distances
r1=sqrt((mu+x)^2);
r2=sqrt((x-(1-mu))^2);

%Compute the Jacobi Energy
C=-(vy^2)/2 + 2*((x^2)/2 + (1-mu)/r1 + mu/r2);
vy = sqrt(2*(2*((x^2)/2 + (1-mu)/r1 + mu/r2) - C2));

sol = solve(C == C2,vy)%,'ReturnConditions',true)