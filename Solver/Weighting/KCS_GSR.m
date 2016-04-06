function [F, result] = KCS_GSR(R, G, Y, opts, image1)
%       F = recSepTVNLM_Wx(R, G, Y, opts, image1)
%       function to recovery SepTV using Split Bregman method with seperable operator 
%       under NLM constrainst
%       Input: 
%           - R, G: sensing matrix size mxn, Phi*X = G'*X*R; or Phi =
%           Kornecker(G', R);
%           - Y: measurement matrix, mxm
%           - opts: parameter
%               + lambda: parameter for gradX and gradY
%               + mu: parameter for ||RFG - Y||_2
%               + nuy: parameter for V-->F
%               + gamma: Parameter for NLM
%               + nInner: maximum loop for Inner Loop, 15-25
%               + nOuter: number of maximum loop for Outer loop, 30-50
%               + tol: tolerance, 6e-7
%   Note: 
%       - The total loop will be: nInner * nOuter, can be adjust between
%         these two parameters. 
%       - NLM is performed in the end of each inner Loop, and was change to
%       boot up the result (0.05-.30 dB), by changing nInner, (When the
%       normal case does not gain), but it increase and then decrease-->
%       take care the "tol" value too.
% Implemented by Thuong Nguyen, 2013, April 03

%% test function
% clear all;
% load testLenna.mat;
%     nInner = 10;    nOuter = 60;       xSize = nInner;
%     tol = 5e-4;   gamma = 2;
%% Assign value;
    load('all_TV_lena_s0.1.mat');
    lambda  = opts.lambda;
    nuy     = opts.nuy ;
    mu      = opts.mu ;
    nInner  = opts.nInner;
    nOuter  = opts.nOuter;
    tol     = opts.tol;
    
%% Initial value    
    [~, N] = size(R);       
    % Reserve memory for the auxillary variables
    Y0 = Y;             F = zeros(N,N);     V = F;			% zeros(N, N);
	x = zeros(N, N);    y = x;              bx = x;
    by = x;             W = x;              Wx = ones(size(F));     Wy = Wx;
 	
	% Decompose the matrix R and G
	[Ur, Lr] = eig(transpose(R)*R);         Ur = Ur';
	[Ug, Lg] = eig(G*transpose(G));         Ug = Ug';
    
	% Decompose differential opertator
    [dx, dy] = gradCal3(F, 3);
    dxt      = dx';                         dyt = dy';
    
	[Ux, Lx] = eig(dxt*dx);                 Ux = Ux';
	[Uy, Ly] = eig(dy*dyt);                 Uy = Uy';
    
	% The denomitor --> recheck
	denoF = nuy + mu*diag(Lr,0)*diag(Lg,0)';
	denoV = nuy + lambda* ( diag(Lx,0)*ones(1,N) + (diag(Ly,0)*ones(1,N))' );
    
    %% Main loop
	errorKPlus = 10;                        outer = 0; 
    
    %% update parameters
    par.R = R;  par.G = G;  
    par.dx=dx;  par.dy=dy;  
    par.Ur=Ur;  par.Lr=Lr;  par.Ug=Ug;  par.Lg=Lg;
    par.Ux=Ux;  par.Lx=Lx;  par.Uy=Ly;  par.Ly=Ly;
    par.init = 0;
    par.V = V;  par.x = x;  par.y = y;  par.F = F;
    par.bx=bx;  par.by=by;  par.W = W;
    par.denoF = denoF;      par.denoV = denoV;    
    nLoop     = nInner;
    
    
    %% main loop - TV -- Initialization 
%     while ((errorKPlus> tol) && outer ~= nOuter)
%         outer = outer +1;       F1 = F;      
%         
%         % Reconstruction
%         % Weight calculation        
%         if(opts.weight_mode ~= 1 && outer > opts.iterWeight)
%             [Wx, Wy] = Weight_Assign(opts.weight_mode, par.F, dx*par.F, par.F*dy, opts.tau1, opts.hWeight);            
% %             figure(); imshow(Wx, []); title(['iter' num2str(outer)]);
%             [res]   = AWTV(Y, lambda, nuy, mu, nLoop, Wx, Wy, par);
%         else
%             [res]   = ATV(Y, lambda, nuy, mu, nLoop, Wx, Wy, par);            
%         end;
%         
%                 
%         F       = res.F;        
%         par.V   = res.V;    par.x   = res.x;    par.y   = res.y;
%         par.bx  = res.bx;   par.by  = res.by;   par.W   = res.W;
%         par.F   = res.F;
%         
%         Y       = Y+Y0-R*F*G;         
% %         % Weight calculation        
% %         if(opts.weightMode ~= 1 && outer > opts.IterEnableWeight)
% %             [Wx, Wy] = Weight_Assign(opts.weightMode, par.F, dx*par.F, par.F*dy, opts.tau1, opts.hWeight);            
% % %             figure(); imshow(Wx, []); title(['iter' num2str(outer)]);
% %         end;
%         errorKPlus = sqrt(sum((F(:) - F1(:)).^2)/sum(F(:).^2)); 
%         
%         if(opts.isShowPSNR == 1)
%             eachErr(outer)  = errorKPlus ;
%             eachPSNR(outer) = csnr(F, image1, 0, 0);
%             eachSSIM(outer) = cal_ssim(F, image1, 0, 0);            
%             display(['   outer' num2str(outer) ', inner' num2str(nInner) ', PSNR:' num2str(eachPSNR(outer)) ...
%                 ', SSIM:' num2str(eachSSIM(outer)) ', errK+:' num2str(eachErr(outer))]);
%         end
% %         figure(1);
% %         imshow(F, []);
%     end
%     
    %% GSR based method 
    load('all_TV_lena_s0.1.mat');
    b = zeros(size(F)); 
    IterNum = 200; 
    Inloop  = 5;
%     Opts = GSR_ParSet(R, F, block_size,image)
    Opts.PatchSize = 8;
    Opts.SlidingDis = 4;
    Opts.Factor = 240;
    Opts.ArrayNo = 60;
    Opts.SearchWin = 20;
    Opts.org = image1;
    if ~isfield(Opts,'IterNum')
    Opts.IterNum = 50;
    end
    if ~isfield(Opts,'mu')
        Opts.mu = 2.5e-3;
    end
    if ~isfield(Opts,'lambda')
        Opts.lambda = 0.082;
    end
    if ~isfield(Opts,'Inloop')
        Opts.Inloop = 200;
    end

    for i = 1:IterNum
        % solving GSR sparse
        F_hat = F; 
        r = F_hat - b; 
        F_bar = GSR_Solver_CS(r, Opts);
        
        % solving subproblem F
        for j = 1:1:Inloop
            % update F
            Yk = mu*(R)'*Y*(G)';
            Fk = Ur*(Yk + nuy*(F_bar + b))*Ug';
            Fk = Fk./denoF;
            F = (Ur)'*Fk*Ug;
            F_hat = F; 
            Y       = Y+Y0-R*F*G;    
            display(['iter out: ' num2str(i) ', in:' num2str(j) ', PSNR ' num2str(csnr(F, image1, 0, 0)) 'dB']);
        end;
        
        b = b - (F_hat - F_bar); 
%         display(['iter ' num2str(i) ', PSNR ' num2str(csnr(F, image1, 0, 0)) 'dB']);
    end
    
    %% results
    result.eacherr = errorKPlus;
    if(opts.isShowPSNR == 1)
        result.eachPSNR = eachPSNR;
        result.eachErr = eachErr;
        result.eachSSIM = eachSSIM;
    end
    
            
     
