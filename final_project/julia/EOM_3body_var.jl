include("EOM_3body.jl")
include("stm_jacobian.jl")

function EOM_3body_var(X,p,t)
	x = X[1];
	y = X[2];
	z = X[3];

	mu = p.mu;

	# from EOM_3body
	Xdot = EOM_3body(X[1:6],p,t);

	Phi = reshape(X[7:42], 6,6);
	# state transition matrix
	A = stm_jacobian(x,y,z,mu);
	Phid = A*Phi;

	state = vcat(Xdot, vec(reshape(Phid, 36, 1)));

	return state

end