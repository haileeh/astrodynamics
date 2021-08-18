using LinearAlgebra

function stm_jacobian(x,y,z,mu)
	# find Jacobian (A) of state transition matrix Phi
	# solve for symmetric matrix of second partial derivatives of U wrt x,y,z
	r1 = sqrt((x+mu)^2 + y^2 + z^2);
	r2 = sqrt((x+mu-1)^2 + y^2 + z^2);

	Omega_xx = (mu - 1)/r1^3 - mu/r2^3 + (0.75*mu*(2*mu + 2*x - 2)^2)/r2^5 - (0.7500*(2*mu + 2*x)^2*(mu - 1))/r1^5 + 1;
	Omega_xy = (1.5*mu*y*(2*mu + 2*x - 2))/r2^5 - (1.5*y*(2*mu + 2*x)*(mu - 1))/r1^5;
	Omega_xz = (1.5*mu*z*(2*mu + 2*x - 2))/r2^5 - (1.5*z*(2*mu + 2*x)*(mu - 1))/r1^5;
	Omega_yx = (1.5*mu*y*(2*mu + 2*x - 2))/r2^5 - (1.5*y*(2*mu + 2*x)*(mu - 1))/r1^5;
	Omega_yy = (mu - 1)/r1^3 - mu/r2^3 + (3*mu*y^2)/r2^5 - (3*y^2*(mu - 1))/r1^5 + 1;
	Omega_yz = (3*mu*y*z)/r2^5 - (3*y*z*(mu - 1))/r1^5;
	Omega_zx = (1.5*mu*z*(2*mu + 2*x - 2))/r2^5 - (1.5*z*(2*mu + 2*x)*(mu - 1))/r1^5;
	Omega_zy = (3*mu*y*z)/r2^5 - (3*y*z*(mu - 1))/r1^5;
	Omega_zz = (mu - 1)/r1^3 - mu/r2^3 + (3*mu*z^2)/r2^5 - (3*z^2*(mu - 1))/r1^5;

	fancyU = [Omega_xx Omega_xy Omega_xz; Omega_yx Omega_yy Omega_yz; Omega_zx Omega_zy Omega_zz]; # matrix
	Delta = zeros(Float64,3,3); #matrix
	Delta[1,2] = 1; Delta[2,1] = -1;
	Ar1 = hcat(zeros(Float64,3,3), Array(1.0I, 3, 3));
	Ar2 = hcat(fancyU, 2.0*Delta);
	A = vcat(Ar1, Ar2);
	return A
end