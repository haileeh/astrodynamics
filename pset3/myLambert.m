function [v1,v2]=myLambert(R1,R2,theta,zinit,T1,mu)

r1 = norm(R1);
r2 = norm(R2);
% solve for A, eqn 5.35 Curtis
A = sind(theta)*sqrt(r1*r2/(1-cosd(theta)));

z(1) = zinit;
i = 1;
while i<10
    %S = @(z) (1/6) - z/120 + z^2/5040 - z^3/362880 + z^4/39916800 - z^5/6227020800;
    %C = @(z) (1/2) - z/24 + z^2/720 - z^3/40320 + z^4/3628800 - z^5/479001600;
    
    if z > 0 % ellipse
        S = @(z) (sqrt(z)-sin(sqrt(z)))/z^(3/2);
        C = @(z) (1-cos(sqrt(z)))/z;
    elseif z < 0
        S = @(z) (sinh(sqrt(-z)) - sqrt(-z))/(-z)^(3/2);
        C = @(z) (cosh(sqrt(-z)) - 1) / -z;
    else % z == 0
        S = @(z) 1/6;
        C = @(z) 1/2;
    end
    
    y = @(z) r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
    F = @(z)(y(z)/C(z))^(3/2)*S(z) + A*sqrt(y(z)) - sqrt(mu)*T1;
    if z(i) ~=0
        dFdt = @(z) (y(z)/C(z))^(3/2)*(1/(2*z)*(C(z)- (3/2)*S(z)/C(z)) + (3/4)*(S(z)^2/C(z))) +...
            A/8*(3*(S(z)/C(z))*sqrt(y(z)) + A*sqrt(C(z)/y(z)));
    else
        dFdt = @(z) sqrt(2)/40*y(0)^(3/2) + A/8*(sqrt(y(0)) + A*sqrt(1/(2*y(0))));
    end
    
    z(i+1) = z(i) - F(z(i))/dFdt(z(i));
    i=i+1;
end

% calculate y
yFinal = y(z(end));  % ND distance

% step 6. lagrange functions
f = 1 - yFinal/r1;
g = A*sqrt(yFinal/mu);
gdot = 1-yFinal/r2;

% step 7. calculate v1 and v2 
v1 = 1/g*(R2-f*R1);
v2 = 1/g*(gdot*R2 - R1);