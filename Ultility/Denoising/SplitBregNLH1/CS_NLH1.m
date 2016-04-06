function [u0,energy,relmse,psnr_n]=CS_NLTV(Fp,Fmask,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   min NLTV(u) s.t. RFu=f, where R is a subsampling matrix, and F is  Fourier transform 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fp: [M,N] complex matrix --------- given Fourier measurement data
%%% For example: Fp=fft2(I).*Fmask/sqrt(M*N), where [M,N]=size(I);
%%% Fmask: [M,N] binary matrix -------------------- Sampling mask, 1:sampled, 0:missed

%%% opts: parameters set. see below for more explanation.
%%% Xiaoqun Zhang: xqzhang@math.ucla.edu

if (nargin<3)
    opts=[];
end

%% Initialization 
% Regularization parameters
if ~isfield(opts,'mu')     opts.mu=10; end %% regularization term scale
if ~isfield(opts,'delta')  opts.delta=1; end %% delta<||A^*A||, since A here is a subsampled Fourier matrix, could be fixed to be 1.
if ~isfield(opts,'nOuter') opts.nOuter=200; end %% Outer Bregman iteration steps.
if ~isfield(opts,'nDenoising') opts.nDenoising=5; end%% 3st level: denoising/regularization level, in general could be fixed to be 10
if ~isfield(opts,'type')  opts.type=1; end %% 1: by bos,2:PBOS. 3:by Linearized Bregman If the without noise (or low), choose 1
if ~isfield(opts,'bTol') opts.bTol=10^-10; end %%stopping criterion on residual std2(Fmask.*fft2(u)/N-Fp), if the noise standard variation is known, can set btol to be sigma
if ~isfield(opts,'xTol') opts.xTol=10^-10; end %%stopping criterion, if the noise standard variation is known, can set btol to be sigma
if ~isfield(opts,'verbose') opts.verbose=0; end %% display message

% Weight parameters 
if ~isfield(opts,'h0')  opts.h0=15; end %% weight filter parater, depends on noise and image standard variation. for example: for barbara [0, 255], h0=20; To be adapted for normalized image
if ~isfield(opts,'nWin') opts.nWin=2;  end    %% patch size [2*nwin+1, 2*nwin+1]
if ~isfield(opts,'nBloc') opts.nBloc=7; end   %% search window size [2*bloc+1, 2*bloc+1]
if ~isfield(opts,'nNeigh') opts.nNeigh=10; end   %% number of best neighbors (real neighbors size: 2*nNeigh+4)
if ~isfield(opts,'nWeightupdate') opts.nWeightupdate=20;  end %0: no weight update, otherwise update steps 
if ~isfield(opts,'denoising_type') opts.denoising_type=1; end %% NLTV denoising algorithm:1: Split bregman, 2: Projection in dual

%% Get start!
[M,N]=size(Fp);
scale=sqrt(M*N);
A=@(x)(fft2(x).*Fmask/scale);
AT=@(x)real(ifft2(x.*Fmask*scale));

u00=AT(Fp);
energy=[];
relmse=[];
psnr_n=[];

u0=u00;
v1=u00;
v0=u00;

tau=1./(opts.delta*opts.mu);
if opts.type==2
    diag=Fmask+1./opts.delta;
end

n=0;
condition=1;
%% Main outer loop: 1st level
while (condition)
  
  if (opts.nWeightupdate) %% update weight
      if(mod(n,opts.nWeightupdate)==0)
        wopts=update_weight(real(v1),opts.h0,opts.nWin,opts.nBloc,opts.nNeigh);          
      end
  end
   
     uold=u0;
     u0_F=A(u0);
     resU=AT(Fp-u0_F); 
    switch opts.type 
        case 1 %% by bregmanized operator splitting
        v0=v0+resU;
        v1= u0+opts.delta*(v0+resU-u00);
        case 2 % pbos 
        v0=v0+resU;
        v1=u0+ real(ifft2(fft2(v0-ifft2(fft2(u0).*Fmask)).*Fmask./diag));        
        case 3   %% by linearized bregman
         v1=v1+opts.delta*resU;       
        case 4 % operator splitting
       v1=u0+ opts.delta* resU;
       otherwise
      disp('Unknown method.')
     end
     
   %  if (opts.denoising_type==1)
 %    u0=denoising_SBNLTV(v1,tau,wopts,opts.nDenoising,tau./10); 
    % u0=denoising_SBNLTV(v1,tau,tau./4,opts.nDenoising,wopts);
     u0=denoising_NLH1(v1,tau,opts.nDenoising,wopts);
    
    % u0=denoising_SBNLTV(v1,tau,tau./10,opts.nDenoising,wopts);
    % else
    % u0=denoising_NLTV_proj(v1,tau,opts.Iters3,wopts);
 %end    
 n=n+1; 
   res=Fp-A(u0);
   energy=[energy;imnorm(res)];    
   relmse= [relmse;sum((u0(:)-uold(:)).^2)./max(sum(u0(:).^2)+sum(uold(:).^2),1)];
     
 if isfield(opts,'I') psnr_n=[psnr_n;psnr(uint8(opts.I),uint8(u0))];end
 if (opts.verbose&mod(n,20)==0) 
     fprintf('\n n=%d,residual std2=%f, relative diff =%f',n-1,energy(n),relmse(n));
   if isfield(opts,'I')  fprintf(',PSNR(u0) =%f',psnr_n(n));end
 end    
 condition=(n<opts.nOuter&& energy(n)>opts.bTol&&relmse(n)>opts.xTol);
end




