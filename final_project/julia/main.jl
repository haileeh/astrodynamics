# main function
using Plots
include("solveLpoints.jl")
include("richardson.jl")
include("differentialControl_zFixed.jl")
include("continuation.jl")
include("plotHaloOrbits.jl")
include("monodromyMatrix.jl")
include("findManifold.jl")

struct Pz
    mu::Float64
end 


plotFlag = 1;
unstableFlag = 0;
# Define constants
mEarth = 5.9722e24; #[kg]
mSun = 1.989e30; #[kg]
mu = mEarth/(mSun+mEarth);
au = 149.5978714*10^6; #km - distance between Sun and Earth

G = (6.67408e-11); # m^3 / kg / s^2
n = -sqrt((mSun+mEarth)*G/((au)*10^3)^3); # rad/s

# amplitude of halo orbit in km
Az = zeros(2,1);
Az[1] = 125000; 
Az[2] = 150000;

# Lagrange points
(L1,L2,L3) = solveLpoints(mu);

# which Lagrange pt?
L1 = convert(Float64,L1);
L2 = convert(Float64,L2);
L=L2;

X0 = zeros(7,6);
t_half = zeros(7,1);
for i in 1:length(Az)
    # Third-order approximation of initial conditions
    x_richardson = richardson(mu,Az[i],n,au,L);
    x_richardson = convert(Array{Float64},x_richardson);

    # Converting from dimensioned coordinates about the Lagrange point
    # to be about the center of mass of the system and ND
    x_3a = [(L*au+x_richardson[1])/au; 0.0; x_richardson[3]/au; 0.0; (L*au-(L*au*n+x_richardson[5])/n)/au; 0.0];
    # fine-tune initial conditions of halo orbit using differential correction
    (X0[i,:],t_half[i])=differentialControl_zFixed(mu,x_3a,3.0); # TODO initial time
end

## continuation!
for k in 1:5
    (X0[2+k,:],t_half[2+k]) = continuation(X0[k,:],X0[k+1,:],mu);
    #println(X0[2+k,:]);
end

## plot halo orbits
if plotFlag==1
    plotHaloOrbits(X0, mu, L1, L2, t_half);
end

## find manifold
p = Pz(mu);
for j in 1:5
    (uD, sD, Y_u, Y_s, X) = monodromyMatrix(X0[j,:],t_half[j],p);

    # find stable manifold
    ic_t = X[1:6,:];
    (state_t,C,tMan)=findManifold(ic_t,p,Y_s,sD); # is the second one corresponding to the stable value?
end
