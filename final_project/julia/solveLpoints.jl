## solveLpoints
# Provided with a value mu, the collinear lagrange points are determined
# L1, L2, L3
using Pkg
Pkg.add("SymPy")
using SymPy
const S = sympy.S
function solveLpoints(mu)

	x = symbols("x");
	y = symbols("y");
	z = symbols("z");

	R1 = sqrt((x+mu)^2 + y^2 + z^2);
	R2 = sqrt((x+mu-1.0)^2 + y^2 + z^2);
	func = 0.5*(x^2 + y^2) + (1.0-mu)/R1 + mu/R2;
	dx = diff(func, x);
	yl = 0.0; zl = 0.0;

	dx_l = dx(y => yl);
	dx_l = dx_l(z => zl);

	sol = solveset(dx_l,(x),domain=S.Reals); 
	sol_elem = elements(sol); # get elements out of finite set
	
	# define points
	tempMax = 0;
	tempMin = 0;
	midVal = 0;
	for i in 1:length(sol_elem)
		if sol_elem[i] > tempMax
			tempMax = sol_elem[i];
		end
		if sol_elem[i] < tempMin
			tempMin = sol_elem[i];
		end
	end
	L2 = tempMax; # find maximum positive value for L2, (around x=1)
	L3 = tempMin; # most negative value for L3, (around x=-1)

	for j in 1:length(sol_elem)
		if sol_elem[j] > tempMin && sol_elem[j] < tempMax
			midVal = sol_elem[j];
		end
	end
	L1 = midVal; # remaining value is L1, (around x=0)


	return L1, L2, L3
end