function u = SB_ITV_fft(g,mu, tol, image)
% Split Bregman Isotropic Total Variation Denoising
%
%   u = arg min_u 1/2||u-g||_2^2 + mu*ITV(u)
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

% load testITVfft

g = g(:);
n = (length(g));
[B, Bt, ~] = DiffOper(sqrt(n));
% mu      = 18; 
lambda  = 0.2;
% prepare for subproblem u
D1 = abs(psf2otf([-1,1],[sqrt(n), sqrt(n)])).^2;
D2 = abs(psf2otf([-1;1],[sqrt(n), sqrt(n)])).^2;
%d1 = lambda*(D1 + D2) + 1;
d1 = lambda*(D1(:) + D2(:)) + 1;
FB0 = fft2(g);
%RHS = speye(n) + lambda * DxyD;
b = zeros(2*n,1);
d = b;
u = g;
err = 1;k = 0;
% tol = 1e-3;

while (err > tol && k < 20)
%     fprintf('it. %g ',k);
    up = u;
    %[u,~] = cgs(speye(n)+BtB, g-lambda*Bt*(b-d),1e-5,100);
    FB  = FB0 + lambda*fft2(Bt*(d - b));
    u   = ifft2(FB./d1);
    u   = real(u);
    
    Bub = B*u+b;
    s = sqrt(Bub(1:n).^2 + Bub(n+1:end).^2);
    d = [max(s-mu/lambda,0).*Bub(1:n)./s ;
        max(s-mu/lambda,0).*Bub(n+1:end)./s ];
    b = Bub-d;
    err = norm(up-u)/norm(u);
    %psnrAll(iter) =  psnr(u(:), image(:), 255); 
%     fprintf('err=%g, PSNR:%g \n',err, psnr(u(:), image(:), 255));
    k = k+1;
end
% imshow(reshape(u,512,512),[]);
% fprintf('Stopped because norm(up-u)/norm(u) <= tol=%.1e\n',tol);
end

function [B, Bt, BtB] = DiffOper(N)
    D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
    D(:,1) = [];
    D(1,1) = 0;
    B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
    Bt = B';
    BtB = Bt*B;
end
