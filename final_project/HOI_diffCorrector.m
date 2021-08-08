function HOI_diffCorrector(tti_state,hoi_state,hoi0_state,mu,tf)
% subscript h identifies velocity components on the nominal path at the HOI
% xdot0, ydot0, zdot0 are the velocity components on the transfer path at
% HOI!
x0 = hoi0_state(1); y0 = hoi0_state(2); z0 = hoi0_state(3);
xdot0 = hoi0_state(4); ydot0 = hoi0_state(5); zdot0 = hoi0_state(6);
xdoth = hoi_state(4); ydoth = hoi_state(5); %zdoth = hoi_state(6);
% x1 = tti_state(1); y1 = tti_state(2); z1 = tti_state(3);
% xdot1 = tti_state(4); ydot1 = tti_state(5); zdot1 = tti_state(6);
% 1 =f
% constants
au = 149.5978714*10^6; %km - distance between Sun and Earth
R = 6371/au; % ND

% call EOM_3body to find xdot1 thru zddot1
p.mu = mu; t = 0;
p.R = R;
% define ha_Star, alphaStar
haStar = 2000/au; %185/au; %ND
alphaStar = -169 * pi/180; %rad

iter = 0;
% HOI equations
DVxy = sqrt( (xdot0 - xdoth)^2 + (ydot0 - ydoth)^2);
beta = atan( ydot0-ydoth / xdot0-xdoth);
%DVz = zd0 - zdh;
while 1
    ic = [x0,y0,z0,xdot0,ydot0,zdot0, reshape(eye(6), 1, 6^2)];
    tspan = linspace(2.5,0,1000); % backwards for stable
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @gammaZero);
    [~,state] = ode45(@EOM_3body_var, tspan, ic, opts, p);
    Phi = reshape(state(end, 7:end), 6, 6);
    x1 = state(end,1);y1=state(end,2);z1=state(end,3);
    xdot1=state(end,4);ydot1=state(end,5);zdot1=state(end,6);
    Xdot1 = EOM_3body(t,state(end,1:6),p);
    
    % TTI equations
    ha = sqrt( (x1-1+mu)^2 + y1^2 + z1^2) - R;
    alpha = atan(y1 / (x1-1+mu) );
    gamma = acos( sqrt( (y1*zdot1 + z1*ydot1)^2 + (z1*xdot1 - (x1-1+mu)*zdot1)^2*( (x1-1+mu)*ydot1 - y1*xdot1)^2) /...
        ( (ha + R)*sqrt(xdot1^2+ydot1^2+zdot1^2)) );
    %lambda = atan( z1 / sqrt((x1-1+mu)^2 + y1^2));
    
    if abs(gamma) < 0.5
        break
    end
    
    % check if simulation can be terminated
    if abs(gamma) < 1e-10
        break;
    end
    
    % expressions for Bij
    dxdot0_dDVxy = ((xdot0 - xdoth)^2 + (ydot0 - ydoth)^2)^(1/2)/(xdot0 - xdoth);
    dydot0_dDVxy = ((xdot0 - xdoth)^2 + (ydot0 - ydoth)^2)^(1/2)/(ydot0 - ydoth);
    dzd0_dDVxy = 0;
    
    dxdot0_dbeta = (xdot0^2*((xdoth - ydot0 + ydoth/xdot0)^2 + 1))/ydoth;
    dydot0_dbeta = (xdoth - ydot0 + ydoth/xdot0)^2 + 1;
    dzd0_dbeta = 0;
    
    dha_dxf = (mu + x1 - 1)/((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2);
    dha_dyf = y1/((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2);
    dha_dzf = z1/((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2);
    
    dalpha_dxf = -y1/((y1^2/(mu + x1 - 1)^2 + 1)*(mu + x1 - 1)^2);
    dalpha_dyf = 1/((y1^2/(mu + x1 - 1)^2 + 1)*(mu + x1 - 1));
    dalpha_dzf = 0;
    
    dgamma_dxf = ((((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*(2*mu + 2*x1 - 2))/(2*((mu + x1 - 1)^2 + y1^2 + z1^2)^(3/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)) - (2*ydot1*(ydot1*(mu + x1 - 1) - xdot1*y1)*(zdot1*(mu + x1 - 1) - xdot1*z1)^2 + 2*zdot1*(ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1))/(2*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)))/(1 - ((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(((mu + x1 - 1)^2 + y1^2 + z1^2)*(xdot1^2 + ydot1^2 + zdot1^2)))^(1/2);
    dgamma_dyf = ((y1*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2))/(((mu + x1 - 1)^2 + y1^2 + z1^2)^(3/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)) - (2*zdot1*(ydot1*z1 + y1*zdot1) - 2*xdot1*(ydot1*(mu + x1 - 1) - xdot1*y1)*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(2*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)))/(1 - ((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(((mu + x1 - 1)^2 + y1^2 + z1^2)*(xdot1^2 + ydot1^2 + zdot1^2)))^(1/2);
    dgamma_dzf = ((z1*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2))/(((mu + x1 - 1)^2 + y1^2 + z1^2)^(3/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)) - (2*ydot1*(ydot1*z1 + y1*zdot1) - 2*xdot1*(ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1))/(2*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)))/(1 - ((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(((mu + x1 - 1)^2 + y1^2 + z1^2)*(xdot1^2 + ydot1^2 + zdot1^2)))^(1/2);
    dgamma_dxdf = ((xdot1*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2))/(((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(3/2)) + (2*y1*(ydot1*(mu + x1 - 1) - xdot1*y1)*(zdot1*(mu + x1 - 1) - xdot1*z1)^2 + 2*z1*(ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1))/(2*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)))/(1 - ((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(((mu + x1 - 1)^2 + y1^2 + z1^2)*(xdot1^2 + ydot1^2 + zdot1^2)))^(1/2);
    dgamma_dydf = ((ydot1*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2))/(((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(3/2)) - (2*z1*(ydot1*z1 + y1*zdot1) + 2*(ydot1*(mu + x1 - 1) - xdot1*y1)*(zdot1*(mu + x1 - 1) - xdot1*z1)^2*(mu + x1 - 1))/(2*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)))/(1 - ((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(((mu + x1 - 1)^2 + y1^2 + z1^2)*(xdot1^2 + ydot1^2 + zdot1^2)))^(1/2);
    dgamma_dzdf = ((zdot1*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2))/(((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(3/2)) - (2*y1*(ydot1*z1 + y1*zdot1) + 2*(ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)*(mu + x1 - 1))/(2*((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)^(1/2)*((mu + x1 - 1)^2 + y1^2 + z1^2)^(1/2)*(xdot1^2 + ydot1^2 + zdot1^2)^(1/2)))/(1 - ((ydot1*z1 + y1*zdot1)^2 + (ydot1*(mu + x1 - 1) - xdot1*y1)^2*(zdot1*(mu + x1 - 1) - xdot1*z1)^2)/(((mu + x1 - 1)^2 + y1^2 + z1^2)*(xdot1^2 + ydot1^2 + zdot1^2)))^(1/2);
    
    % partial derivatives
    dha_dDVxy = [dha_dxf, dha_dyf, dha_dzf]*Phi(1:3,4:6)*[dxdot0_dDVxy;dydot0_dDVxy;dzd0_dDVxy];
    dgamma_dDVxy = [dgamma_dxf,dgamma_dyf,dgamma_dzf,dgamma_dxdf,dgamma_dydf,dgamma_dzdf]*Phi(1:6,4:6)*[dxdot0_dDVxy;dydot0_dDVxy;dzd0_dDVxy];
    dha_dsig1 = [dha_dxf, dha_dyf, dha_dzf]*[xdot1;ydot1;zdot1];
    dgamma_dsig1 = [dgamma_dxf,dgamma_dyf,dgamma_dzf,dgamma_dxdf,dgamma_dydf,dgamma_dzdf]*Xdot1;
    dha_dbeta = [dha_dxf, dha_dyf, dha_dzf]*Phi(1:3,4:6)*[dxdot0_dbeta; dydot0_dbeta; dzd0_dbeta];
    dgamma_dbeta = [dgamma_dxf,dgamma_dyf,dgamma_dzf,dgamma_dxdf,dgamma_dydf,dgamma_dzdf]*Phi(1:6,4:6)*[dxdot0_dbeta; dydot0_dbeta; dzd0_dbeta];
    dalpha_dDVxy = [dalpha_dxf, dalpha_dyf, dalpha_dzf]*Phi(1:3,4:6)*[dxdot0_dDVxy;dydot0_dDVxy;dzd0_dDVxy];
    dalpha_dbeta = [dalpha_dxf, dalpha_dyf, dalpha_dzf]*Phi(1:3,4:6)*[dxdot0_dbeta; dydot0_dbeta; dzd0_dbeta];
    % matrices
    B11 = dha_dDVxy - (dha_dsig1 / dgamma_dsig1) * dgamma_dDVxy;
    B12 = dha_dbeta - (dha_dsig1 / dgamma_dsig1) * dgamma_dbeta;
    B21 = dalpha_dDVxy - (dha_dsig1 / dgamma_dsig1) * dgamma_dDVxy;
    B22 = dalpha_dbeta - (dha_dsig1 / dgamma_dsig1) * dgamma_dbeta;
    
    B = [B11, B12; B21, B22];
    deltaState = B\[haStar - ha; alphaStar - alpha]; %[del Delta Vxy; del beta]
    % add delta state to previous values
    DVxy = DVxy + deltaState(1); % changes xdot0, ydot0
    beta = beta + deltaState(2); % changes xdot0, ydot0
    
    ydot0 = ydoth + sqrt( DVxy^2*tan(beta)^2/(1+tan(beta)^2));
    xdot0 = xdoth + (ydot0-ydoth)/tan(beta);
    
    iter = iter +1;
    if iter > 20
        break
    end
    
end
% end
%
% function zdot = rhs(t,z,p)
% %X0 = {x0,y0,z0,Delta Vxy, beta, Delta Vz} HOI
% %X1 = {h_a,alpha,lambda,gamma,V1,zeta} TTI
% end

% syms xdot0 xdoth ydot0 ydoth zd0 zdh real
% syms x1 y1 z1 mu real
% syms xdot1 ydot1 zdot1 R real

% equations
% DVxy = sqrt( (xdot0 - xdoth)^2 + (ydot0 - ydoth)^2);
% dDVxy_dxdot0 = simplify(diff(DVxy,xdot0)); dxdot0_dDVxy = 1 / dDVxy_dxdot0;
% dDVxy_dydot0 = simplify(diff(DVxy,ydot0)); dydot0_dDVxy = 1 / dDVxy_dydot0;
% dDVxy_dzd0 = simplify(diff(DVxy,zd0)); dzd0_dDVxy = 0;

% beta = atan( ydot0-ydoth / xdot0-xdoth);
% dbeta_dxdot0 = simplify(diff(beta,xdot0)); dxdot0_dbeta = 1 / dbeta_dxdot0;
% dbeta_dydot0 = simplify(diff(beta,ydot0)); dydot0_dbeta = 1 / dbeta_dydot0;
% dbeta_dzd0 = simplify(diff(beta,zd0)); dzd0_dbeta = 0;

%%%DVz = zd0 - zdh;
% ha = sqrt( (x1-1+mu)^2 + y1^2 + z1^2) - R;
% dha_dxf = simplify(diff(ha,x1));
% dha_dyf = simplify(diff(ha,y1));
% dha_dzf = simplify(diff(ha,z1));

% alpha = atan(y1 / (x1-1+mu) );
% dalpha_dxf = simplify(diff(alpha,x1));
% dalpha_dyf = simplify(diff(alpha,y1));
% dalpha_dzf = simplify(diff(alpha,z1));

%%%lambda = atan( z1 / sqrt((x1-1+mu)^2 + y1^2));
% gamma = acos( sqrt( (y1*zdot1 + z1*ydot1)^2 + (z1*xdot1 - (x1-1+mu)*zdot1)^2*( (x1-1+mu)*ydot1 - y1*xdot1)^2) /...
%     ( (ha + R)*sqrt(xdot1^2+ydot1^2+zdot1^2)) );
% dgamma_dxf = simplify(diff(gamma,x1));
% dgamma_dyf = simplify(diff(gamma,y1));
% dgamma_dzf = simplify(diff(gamma,z1));
% dgamma_dxdf = simplify(diff(gamma,xdot1));
% dgamma_dydf = simplify(diff(gamma,ydot1));
% dgamma_dzdf = simplify(diff(gamma,zdot1));