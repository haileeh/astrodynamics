function jacobiConst(x,v,mu)

	# computes Jacobi Energy for a given state (x,v)

	#the distances
	r1=sqrt( (mu+x[1])^2 + (x[2])^2 + (x[3])^2 );
	r2=sqrt( (x[1]-(1-mu))^2 + (x[2])^2 + (x[3])^2 );

	#Compute the Jacobi Energy
	C = -(v[1]^2 + v[2]^2 + v[3]^2) + 2*( (1-mu)/r1 + mu/r2) + x[1]^2 + x[1]^2;

	return C
end
