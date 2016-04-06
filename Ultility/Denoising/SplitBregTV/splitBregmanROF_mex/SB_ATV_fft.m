function u = SB_ATV_fft(g, mu, tol, image)
% Split Bregman Anisotropic Total Variation Denoising
%
%   u = arg min_u 1/2||u-g||_2^2 + mu*ATV(u)
%   
%   g : noisy image
%   mu: regularisation parameter
%   u : denoised image
%   tol: stoppin eror
%
% Refs:
%  *Goldstein and Osher, The split Bregman method for L1 regularized problems
%   SIAM Journal on Imaging Sciences 2(2) 2009
%  *Micchelli et al, Proximity algorithms for image models: denoising
%   Inverse Problems 27(4) 2011
%
% Benjamin Trémoulhéac
% University College London
% b.tremoulheac@cs.ucl.ac.uk
% April 2012
g = g(:);
n = length(g);
[B, Bt, BtB] = DiffOper(sqrt(n));
lambda  = 0.2;
% prepare for subproblem u
D1 = abs(psf2otf([-1,1],[sqrt(n), sqrt(n)])).^2;
D2 = abs(psf2otf([-1;1],[sqrt(n), sqrt(n)])).^2;
%d1 = lambda*(D1 + D2) + 1;
d1  = lambda*(D1(:) + D2(:)) + 1;
FB0 = fft2(g);

b = zeros(2*n,1);
d = b;
u = g;
err = 1;k = 1;
% tol = 1e-3;
% lambda = 1;
k = 0; 
while (err > tol && k < 20)
%     fprintf('it. %g ',k);    
    % u-subproblem by gradient descent
    up = u;
   % [u,~] = cgs(speye(n)+BtB, g-lambda*Bt*(b-d),1e-5,100); 
    FB  = FB0 + lambda*fft2(Bt*(d - b));
    u   = ifft2(FB./d1);
    u   = real(u);
    % d-subproblme by shirnkage 
    Bub = B*u+b;
    d = max(abs(Bub)-mu/lambda,0).*sign(Bub);
    
    % update different
    b = Bub-d;
    err = norm(up-u)/norm(u);
%     fprintf('err=%g, PSNR:%g \n',err, psnr(u, image, 255));
    k = k+1;    
end
% fprintf('Stopped because norm(up-u)/norm(u) <= tol=%.1e\n',tol);
end

function [B, Bt, BtB] = DiffOper(N)

% function to calculate different paprameter

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;
end
