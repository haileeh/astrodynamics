function C=jacobiConst_2D(y,v,mu)

%File computes Jacobi Energy for a given state (x,v)

%the distances
r1=sqrt((mu+y(1,1))^2);
r2=sqrt((y(1,1)-(1-mu))^2);

%Compute the Jacobi Energy
C=-(v(2,1)^2)/2 + 2*((y(1,1)^2)/2 + (1-mu)/r1 + mu/r2);

