using LinearAlgebra
include("EOM_3body_var.jl")

struct Ptcph
    mu::Float64
    ic::Array{Float64} 
end 

function condition(z,t,integrator)
    z[2]
end

function affect!(integrator)
    terminate!(integrator)
end

function monodromyMatrix(X0,t_half,p)
	# assume single IC passed in

	ic = vcat(X0, vec(reshape(Array(1.0I, 6, 6),36,1)));

    tspan = (0.0, 2.0*t_half);
    #[T,X] = ode45(@EOM_3body_var, tspan, ic, opts2, p);
    cb = ContinuousCallback(condition,affect!)
    prob = ODEProblem(EOM_3body_var,ic,tspan,p);
    sol = solve(prob, Tsit5(), reltol=1e-8,abstol=1e-8,callback=cb);
    T = sol.t;
    X = sol;

    #lastXi = lastindex(X,2); 
    Phi = zeros(6,6,length(T));
    for i in 1:length(T)
        Phi[:,:,i] = reshape(X[7:42,i], 6, 6);
    end
    Phif = Phi[:,:,length(T)];
    
    V = eigvecs(Phif);
    D = eigvals(Phif);
    Y_u = V[:,1]; #unstable
    Y_s = V[:,2]; #stable
end