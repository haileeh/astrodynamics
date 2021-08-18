function richardson(mu,Az,n,au,L1)

    # Richardson's third-order analytical solution 

    # constants
    n_13 = 1; # bifurcation; values are 1 or 3

    # distance of secondary body from the Lagrange point
    rL = au*abs(1-L1); #km 

    Az = Az/rL; # normalize with respect to rL
    gammaL = rL /au;

    # richardson 1980 analytic construction paper,eqn 8a
    if L1 < 1 #L1
        c2 = (1/gammaL^3)* (mu + (1-mu)*gammaL^3/(1-gammaL)^3);
        c3 = (1/gammaL^3)* (mu - (1-mu)*gammaL^4/(1-gammaL)^4);
        c4 = (1/gammaL^3)* (mu + (1-mu)*gammaL^5/(1-gammaL)^5);
    else #L2
        c2 = (1/(gammaL^3))*( (-1)^2*mu + (-1)^2 * (((1-mu)*gammaL^3)/((1+gammaL)^3)));
        c3 = (1/(gammaL^3))*( (-1)^3*mu + (-1)^3 * (((1-mu)*gammaL^4)/((1+gammaL)^4)));
        c4 = (1/(gammaL^3))*( (-1)^4*mu + (-1)^4 * (((1-mu)*gammaL^5)/((1+gammaL)^5)));
    end
    # linearized frequency lambda by solving
    lambda2 = (-(c2-2) + sqrt((c2-2)^2 + 4*1*(c2-1)*(1+2*c2)))/(2*1);
    lambda = sqrt(lambda2);
    k = 2*lambda / (lambda^2 + 1 - c2);

    d1 = 3*lambda^2 / k * (k*(6*lambda^2 -1) - 2*lambda);
    d2 = 8*lambda^2 / k * (k*(11*lambda^2 -1) - 2*lambda);

    a21 = 3*c3*(k^2-2)/ (4*(1+2*c2));
    a22 = 3*c3 / (4*(1+2*c2));
    a23 = -3*c3*lambda/(4*k*d1) * (3*k^3*lambda - 6*k*(k-lambda)+4);
    a24 = -3*c3*lambda/(4*k*d1) * (2+3*k*lambda);

    d21 = -c3/(2*lambda^2);
    d31 = 3/(64*lambda^2)*(4*c3*a24+c4);
    d32 = 3/(64*lambda^2) * (4*c3*(a23-d21) + c4*(4+k^2));

    b21 = -3*c3*lambda/(2*d1)*(3*k*lambda-4);
    b22 = 3*c3*lambda / d1;
    b31 = 3/(8*d2)* (8*lambda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2)) + (9*lambda^2+1+2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2)));
    b32 = 1/d2* (9*lambda*(c3*(k*b22 + d21 - 2*a24) - c4) + (3/8)*(9*lambda^2 + 1 + 2*c2)*(4*c3*(k*a24-b22)+k*c4));

    a31 = -9*lambda/(4*d2) * (4*c3*(k*a23-b21) + k*c4*(4+k^2)) + (9*lambda^2+1-c2)/(2*d2)*(3*c3*(2*a23-k*b21)+c4*(2+3*k^2));
    a32 = -1/d2*( 9*lambda/4 * (4*c3*(k*a24-b22) + k*c4) + 3/2*(9*lambda^2+1-c2)*(c3*(k*b22+d21-2*a24)-c4));
    # frequency corrections s1 and s2
    s1 = 1/(2*lambda*(lambda*(1+k^2)-2*k)) * (3/2*c3*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21) - 3/8*c4*(3*k^4-8*k^2+8));
    s2 = 1/(2*lambda*(lambda*(1+k^2)-2*k)) * (3/2*c3*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+5*d21) + 3/8*c4*(12-k^2));

    a1 = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-k^2);
    a2 = 3/2*c3*(a24-2*a22) + 9/8*c4;

    l1 = a1+2*lambda^2*s1;
    l2 = a2+2*lambda^2*s2;
    Delta = lambda^2-c2; # nomalized units
    #Azt = sqrt((-Delta - l1*(Ax_min*rL).^2)./l2);

    if L1<1 #L1
        Ax = -sqrt((-l2*Az^2-Delta)/l1);
        Ay = k*Ax;
    else #L2
        Ax_min = sqrt(abs(Delta/l1));
        Ax = -Ax_min;
        Ay = k*Ax;
    end

    omega = 1 + s1*Ax^2 + s2*Az^2; # of the correct order O(A_x^n)
    # solve
    t = 0; 
    s = n*t;
    phi = 0; #determines the family of orbits
    psi = phi + n_13*pi/2; 
    tau = omega*s;
    tau1 = lambda*tau + phi; #n=[rad/s], t =s, omega = [km^2]

    x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Ax^2)*cos(2*tau1) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);
    y = k*Ax*sin(tau1) + (b21*Ax^2-b22*Az^2)*sin(2*tau1)+ (b31*Ax^3-b32*Ax*Az^2)*sin(3*tau1);
    deln = 2 - n_13;
    z = deln*(Az*cos(tau1) + d21*Ax*Az*(cos(2*tau1) - 3) + (d32*Az*Ax^2 -d31*Az^3)*cos(3*tau1));

    # initial conditions for a halo orbit in an Lagrange point centered coordinate frame
    # normalized by the Lagrange point secondary distance (rL)
    ydot0 = k*Ax*(lambda*cos(tau1)) + (b21*Ax^2 - b22*Az^2)*(2*lambda*cos(2*tau1)) + (b31*Ax^3 - b32*Ax*Az^2)*(3*lambda*cos(3*tau1));

    # The coordinates are with respect to the Lagrange point as origin. 
    # rL is the distance between the secondary body and the Lagrange point. The
    # conversion below is converting it from ND to km.
    rE = au-L1*au; 
    X0 = [x*rE 0 z*rE 0 ydot0*n*rE 0];

    return X0

end 