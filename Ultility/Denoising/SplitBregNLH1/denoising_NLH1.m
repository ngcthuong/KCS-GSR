function u1=denoising_NLH1(u0,mu,inner,wopts)
 [Ny,Nx]=size(u0);
 %mu=1./mu; % to be consistent with the code: mu is on the fidelity term in the denosing code
display_messages = 0;

 VecGeneralParameters = [display_messages; Ny; Nx; wopts.m; wopts.w; wopts.NbNeigh; mu; inner;];
u1 = NLH1_mex(single(u0),single(wopts.W),int32(wopts.Y),single(VecGeneralParameters));
