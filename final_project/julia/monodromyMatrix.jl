using LinearAlgebra
include("EOM_3body_var.jl")

struct Ptcph
    mu::Float64
    ic::Array{Float64} 
end 

function monodromyMatrix(X0,t_half,p)
	# assume single IC passed in

	ic = vcat(X0, vec(reshape(Array(1.0I, 6, 6),36,1)));
    tspan = (0.0, 2.0*t_half);

    prob = ODEProblem(EOM_3body_var,ic,tspan,p);
    sol = solve(prob, Tsit5(), reltol=1e-16,abstol=1e-16,adaptive=false,dt=0.0315); # force the time steps, want 100
    T = sol.t;
    X = sol;

    #lastXi = lastindex(X,2);
    Phi = zeros(6,6,length(T));
    for i in 1:length(T)
        Phi[:,:,i] = reshape(X[7:42,i], 6, 6);
    end
    Phif = Phi[:,:,length(T)];
    
    V = eigvecs(Phif); # columns associated with eigvals
    D = eigvals(Phif);
    # unstable vs stable depends on eigenvalue
    # D element > 1 is unstable, D element < 1 is stable
    uIdx = 0;
    sIdx = 0;
    uD = 0.0;
    sD = 100.0;
    for i in 1:length(D)
        if imag(D[i]) == 0.0
            if real(D[i]) > uD && real(D[i]) > 1.0
                uD = real(D[i]);
                uIdx = i;
            elseif real(D[i]) < sD 
                sD = real(D[i]);
                sIdx = i;
            end
        end
    end

    Y_u = real(V[:,uIdx]); #unstable
    Y_s = real(V[:,sIdx]); #stable

    return uD,sD, Y_u, Y_s, X
end