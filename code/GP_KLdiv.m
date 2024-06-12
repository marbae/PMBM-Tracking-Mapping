function [KLdiv] = GP_KLdiv(GP1,GP2)

% Function that computes the KL-divergence between two
% Gamma-Gaussian-inverse Wishart (GP) distributions
%
% p(g,x,X) = Gam(g;a1,b1)N(x;m1,P1)IW(X;v1,V1)
% q(g,x,X) = Gam(g;a2,b2)N(x;m2,P2)IW(X;v2,V2)
%
% D_KL = D(p||q) = int p ln(p/q) dx

% Dimension of kinematical state
n_x = 6;
cov1 = GP1.cov;
cov2 = GP2.cov;
m1 = GP1.state;
m2 = GP2.state;

m1(3) = [];
cov1(:,3) = []; cov1(3,:) = [];
m2(3) = [];
cov2(:,3) = []; cov2(3,:) = [];

KLdiv = (...
    -0.5*log(det(cov1))+0.5*log(det(cov2))...
    -0.5*n_x+0.5*(m1-m2)'*(cov2\(m1-m2))...
    +0.5*trace(cov2\cov1)...
    );

KLdiv_G = GP1.alpha*log(GP1.beta)-GP2.alpha*log(GP2.beta)+gammaln(GP2.alpha)-gammaln(GP1.alpha)...
    +(GP1.alpha-GP2.alpha)*(psi(0,GP1.alpha)-log(GP1.beta))+GP1.alpha*(GP2.beta/GP1.beta-1);

KLdiv = KLdiv+KLdiv_G;

end