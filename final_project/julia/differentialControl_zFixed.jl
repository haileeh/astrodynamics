using DifferentialEquations
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

function differentialControl_zFixed(mu,X0,t_half)
    ## differential correction 
    # Given a guess for the initial conditions, differential correction finds
    # the correct initial conditions that will produce a trajectory to the
    # desired final state

    # For halo orbit, I set the desired final state as the initial state
    # Use symmetry about the x-z plane and periodicity

    # Yields initial conditions and t_half (or just t)
    # example that works

    # X0 be a vector, change this in main

    global it = 0;
    tol = 1e-8;
    tf = t_half;
    while it <= 10
        tspan = (0,tf); #range(0,stop=tf,length=2000);
        Phi0 = Array(1.0I, 6, 6);
        ic = vcat(X0, vec(reshape(Phi0, 36, 1)));
        p = Ptcph(mu,ic);

        cb = ContinuousCallback(condition,affect!)
        prob = ODEProblem(EOM_3body_var,ic,tspan,p);
        sol = solve(prob, Tsit5(), reltol=1e-8,abstol=1e-8,callback=cb);
        T = sol.t;
        X = sol;
        #opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events', @y_cross); 
        #[T,X] = ode45(@EOM_3body_var, tspan, ic, opts, p);

        lastXi = lastindex(X,2); # correct dimension
        x_t1 = X[1:6,lastXi];
        del_x1 = X0 - x_t1; # desired is initial state
        
        global it = it + 1;
        if it > 10 # done yet?
            println("INCOMPLETE")
            break
        end
        
        # check tolerance
        if abs(x_t1[4]) < tol && abs(x_t1[6]) < tol
            println("converged")
            t_half = T[lastXi];
            break;
        else
            # continue
            Phi = reshape(X[7:42,lastXi], 6, 6);
            Xdot = EOM_3body_var(X[1:42,lastXi],p,0.0);           

            # Sanchez
            rhsMat = [Phi[2,1] Phi[2,5] Xdot[2]; Phi[4,1] Phi[4,5] Xdot[4]; Phi[6,1] Phi[6,5] Xdot[6]];
            if abs(det(inv(rhsMat))) < 1e-16
                keyboard
            end
            deltaState = rhsMat \ [0; del_x1[4]; del_x1[6]];
            X0[1] = X0[1] + deltaState[1]; #x0
            X0[5] = X0[5] + deltaState[2]; #ydot0
            tf = tf+deltaState[3];
        end
    end
    return X0, t_half

end