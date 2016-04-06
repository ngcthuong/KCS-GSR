function u1=denoising_SBNLTV(u0,mu,lambda,inner,wopts)
 [Ny,Nx]=size(u0);
 %mu=1./mu; % to be consistent with the code: mu is on the fidelity term in the denosing code
display_messages = 0;

 VecGeneralParameters = [ display_messages; Ny; Nx; wopts.m; wopts.w; wopts.NbNeigh; lambda; mu; inner;];

d = zeros(Nx,Ny,wopts.NbNeigh);
b = zeros(Nx,Ny,wopts.NbNeigh);

u1 = SBNLTV_mex(single(u0),single(d),single(b),single(u0),...
    single(wopts.W),int32(wopts.Y),single(VecGeneralParameters));
