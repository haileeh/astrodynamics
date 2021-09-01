include("jacobiConst.jl")

function findManifold(ic_t,p,Y_s,eigD)
    # ic_t of size 6 x n
    n = size(ic_t,2);

    if Y_s[1] > 0.0
        #neg eps
        eps = -0.0001;
    else
        # pos eps
        eps = 0.0001;
    end

    tf = 2.5;
    #tf = 10;
    if eigD < 1.0
        # use reverse direction
        tspan = (tf,0.0);#linspace(tf,0,1000);
    else
        tspan = (0.0,tf);#linspace(0,tf,1000);
    end
    # determine dt for 1000 time steps
    dtFixed = (tf-0.0)/998.5;
    state_t = zeros(1000,6,100);
    tMan = zeros(1000,100);

    for i in 1:n
       ic = ic_t[:,i] + eps*Y_s;
       #[tMan(:,i),state_t(:,:,i)] = ode45(@EOM_3body, tspan, ic, opts, p);
       prob = ODEProblem(EOM_3body,ic,tspan,p);
       sol = solve(prob, Tsit5(), reltol=1e-12,abstol=1e-12,adaptive=false,dt=dtFixed);
       tMan[:,i] = sol.t;
       state_t[:,:,i] = sol';
    end
    C = zeros(n,1);

    for k in 1:n
        C(k)=jacobiConst(state_t[1,1:3,k],state_t[1,4:6,k],p.mu);
    end

    return state_t,C,tMan

end