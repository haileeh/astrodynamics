%function HOI_diffCorrector()
syms xd0 xdh yd0 ydh zd0 zdh real
syms xf yf zf mu real
syms xdf ydf zdf R real

% equations
DVxy = sqrt( (xd0 - xdh)^2 + (yd0 - ydh)^2);
dDVxy_dxd0 = simplify(diff(DVxy,xd0)); dxd0_dDVxy = 1 / dDVxy_dxd0;
dDVxy_dyd0 = simplify(diff(DVxy,yd0)); dyd0_dDVxy = 1 / dDVxy_dyd0;
dDVxy_dzd0 = simplify(diff(DVxy,zd0)); dzd0_dDVxy = 0;

beta = atan( yd0-ydh / xd0-xdh);
dbeta_dxd0 = simplify(diff(beta,xd0)); dxd0_dbeta = 1 / dbeta_dxd0;
dbeta_dyd0 = simplify(diff(beta,yd0)); dyd0_dbeta = 1 / dbeta_dyd0; 
dbeta_dzd0 = simplify(diff(beta,zd0)); dzd0_dbeta = 0;

%%%DVz = zd0 - zdh;
ha = sqrt( (xf-1+mu)^2 + yf^2 + zf^2) - R;
dha_dxf = simplify(diff(ha,xf));
dha_dyf = simplify(diff(ha,yf));
dha_dzf = simplify(diff(ha,zf));

alpha = atan(yf / (xf-1+mu) );
dalpha_dxf = simplify(diff(alpha,xf));
dalpha_dyf = simplify(diff(alpha,yf));
dalpha_dzf = simplify(diff(alpha,zf));

%%%lambda = atan( zf / sqrt((xf-1+mu)^2 + yf^2));
gamma = acos( sqrt( (yf*zdf + zf*ydf)^2 + (zf*xdf - (xf-1+mu)*zdf)^2*( (xf-1+mu)*ydf - yf*xdf)^2) /...
    ( (ha + R)*sqrt(xdf^2+ydf^2+zdf^2)) );
dgamma_dxf = simplify(diff(gamma,xf));
dgamma_dyf = simplify(diff(gamma,yf));
dgamma_dzf = simplify(diff(gamma,zf));
dgamma_dxdf = simplify(diff(gamma,xdf));
dgamma_dydf = simplify(diff(gamma,ydf));
dgamma_dzdf = simplify(diff(gamma,zdf));

% partial derivatives
dha_dDVxy = [dha_dxf, dha_dyf, dha_dzf]*Phi(1:3,4:6)*[dxd0_dDVxy;dyd0_dDVxy;dzd0_dDVxy];
dgamma_dDVxy = [dgamma_dxf,dgamma_dyf,dgamma_dzf,dgamma_dxdf,dgamma_dydf,dgamma_dzdf]*Phi(1:6,4:6)*[dxd0_dDVxy;dyd0_dDVxy;dzd0_dDVxy];
dha_dsig1 = [dha_dxf, dha_dyf, dha_dzf]*[xdf;ydf;zdf];
dgamma_dsig1 = [dgamma_dxf,dgamma_dyf,dgamma_dzf,dgamma_dxdf,dgamma_dydf,dgamma_dzdf]*[xdf;ydf;zdf;xddf;yddf;zddf];
dha_dbeta = [dha_dxf, dha_dyf, dha_dzf]*Phi(1:3,4:6)*[dxd0_dbeta; dyd0_dbeta; dzd0_dbeta];
dgamma_dbeta = [dgamma_dxf,dgamma_dyf,dgamma_dzf,dgamma_dxdf,dgamma_dydf,dgamma_dzdf]*Phi(1:6,4:6)*[dxd0_dbeta; dyd0_dbeta; dzd0_dbeta];
dalpha_dDVxy = [dalpha_dxf, dalpha_dyf, dalpha_dzf]*Phi(1:3,4:6)*[dxd0_dDVxy;dyd0_dDVxy;dzd0_dDVxy];
dalpha_dbeta = [dalpha_dxf, dalpha_dyf, dalpha_dzf]*Phi(1:3,4:6)*[dxd0_dbeta; dyd0_dbeta; dzd0_dbeta];
% matrices
B11 = dha_dDVxy - (dha_dsig1 / dgamma_dsig1) * dgamma_dDVxy;
B12 = dha_dbeta - (dha_dsig1 / dgamma_dsig1) * dgamma_dbeta;
B21 = dalpha_dDVxy - (dha_dsig1 / dgamma_dsig1) * dgamma_dDVxy;
B22 = dalpha_dbeta - (dha_dsig1 / dgamma_dsig1) * dgamma_dbeta;