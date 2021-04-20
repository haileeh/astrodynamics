%--------------------------------------------------------------------------
%
% Stumpff: computes values for the Stumpff functions C1, C2 and C3
%
% Input:
%   E2       Square of eccentric anomaly (E2=E*E) in [rad^2]
%
% Outputs:
%   c1       Value of C1 = sin(E)/E
%   c2       Value of C2 = (1-cos(E))/(E*E)
%   c3       Value of C3 = (E-sin(E))/(E^3)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [c1, c2, c3] = Stumpff (E2)

epsilon = 100*eps;

c1 = 0;
c2 = c1;
c3 = c2; 

add = 1;
n = add;

while(1)
    c1 = c1 + add;
    add = add/(2*n);
    c2 = c2 + add;
    add = add/(2*n+1);
    c3 = c3 + add;
    add = add*(-E2); 
    n = n + 1;    
    if (abs(add) < epsilon)
        break
    end
end
