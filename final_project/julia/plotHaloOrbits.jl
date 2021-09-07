using Plots
gr()
struct Py
    mu::Float64
end 

function plotHaloOrbits(X0, mu, L1, L2, t_half)
	plot([1.0-mu],[0.0],[0.0], linecolor=:blue, markershape=:circle, lw=2, lab="Earth",xlabel="X",ylabel="Y",zlabel="Z") # Earth
    plot!([L1],[0.0],[0.0],linecolor=:purple, markercolor=:purple, markershape=:hexagon, lw=2, lab="L1") # L1
    plot!([L2],[0.0],[0.0],linecolor=:purple, markercolor=:purple, markershape=:square,lw=2, lab="L2") # L2
    for j in 1:size(X0,1)
        ic = X0[j,1:6];
        p = Py(mu);
        tspan = (0.0,2.0*t_half[j]);

        prob = ODEProblem(EOM_3body,ic,tspan,p);
        sol = solve(prob, Tsit5(), reltol=1e-8,abstol=1e-8,adaptive=false,dt=0.03);
        
        if j==1
        	plot!(sol[1,:],sol[2,:],sol[3,:],linecolor=:black,lab="Halo Orbits")
        else
        	plot!(sol[1,:],sol[2,:],sol[3,:],linecolor=:black,lab="")
        end
    end
    plot!(legend=:outertopright,show=true)
    title!("Halo orbit family about L2")
    savefig("/Users/hettrick/Documents/GitHub/astrodynamics/final_project/julia/halo.svg")
    # rotates about lagrange point from y-z projection
end

function plotStableManifolds(X, state_t, mu, L1, L2,h)
    plot([1.0-mu],[0.0],[0.0], linecolor=:blue, markershape=:circle, lw=2, lab="Earth",xlabel="X",ylabel="Y",zlabel="Z"); # Earth
    plot!([L1],[0.0],[0.0],linecolor=:purple, markercolor=:purple, markershape=:hexagon, lw=2, lab="L1"); # L1
    plot!([L2],[0.0],[0.0],linecolor=:purple, markercolor=:purple, markershape=:square,lw=2, lab="L2"); # L2
    plot!(X[1,:],X[2,:],X[3,:],linecolor=:black,lw=3,lab="Halo Orbit");
 
    for j in 1:3:size(state_t,3) 
        if j==1
            plot!(state_t[:,1,j],state_t[:,2,j],state_t[:,3,j],linecolor=:green,lab="manifold"); # add label for first one
        else
            plot!(state_t[:,1,j],state_t[:,2,j],state_t[:,3,j],linecolor=:green,lab="");
        end
    end
    plot!(xlims = (0.985,1.035), ylims=(-0.0055,0.0045), zlims=(-0.005,0.005) )
    plot!(legend=:outertopright,show=true);
    title!("Stable Manifolds");
    savefig("/Users/hettrick/Documents/GitHub/astrodynamics/final_project/julia/manifold$h.svg")
    # rotates about lagrange point from y-z projection
end

function plotUnstableManifolds(X, state_t, mu, L1, L2,h)
    plot([1.0-mu],[0.0],[0.0], linecolor=:blue, markershape=:circle, lw=2, lab="Earth",xlabel="X",ylabel="Y",zlabel="Z"); # Earth
    plot!([L1],[0.0],[0.0],linecolor=:purple, markercolor=:purple, markershape=:hexagon, lw=2, lab="L1"); # L1
    plot!([L2],[0.0],[0.0],linecolor=:purple, markercolor=:purple, markershape=:square,lw=2, lab="L2"); # L2
    plot!(X[1,:],X[2,:],X[3,:],linecolor=:black,lw=3,lab="Halo Orbit");
 
    for j in 1:3:size(state_t,3) 
        if j==1
            plot!(state_t[:,1,j],state_t[:,2,j],state_t[:,3,j],linecolor=:red,lab="manifold"); # add label for first one
        else
            plot!(state_t[:,1,j],state_t[:,2,j],state_t[:,3,j],linecolor=:red,lab="");
        end
    end
    plot!(xlims = (0.985,1.035), ylims=(-0.0055,0.0045), zlims=(-0.005,0.005) )
    plot!(legend=:outertopright,show=true);
    title!("Unstable Manifolds");
    savefig("/Users/hettrick/Documents/GitHub/astrodynamics/final_project/julia/manifold$h.svg")
    # rotates about lagrange point from y-z projection
end