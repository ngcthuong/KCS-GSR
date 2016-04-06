function [F result] = recSepWTVNL(R, G, Y, opts, image1, mode)
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
    lambda      = opts.lambda;           nuy    = opts.nuy ;
    mu          = opts.mu ;              nInner = opts.nInner;
    nOuter      = opts.nOuter;           tol    = opts.tol;
    par.gamma   = opts.gamma;
    
    Options.kernelratio = 6 ;%: %Radius of local Patch (default 3)
    Options.windowratio = 3;%: Radius of neighbourhood search window (default 3)
    Options.filterstrength = opts.hReg;
    
%% Initial value    
    tic
    [M, N] = size(R);       
    % Reserve memory for the auxillary variables
    Y0      = Y;             F = zeros(N,N);     V   = F;			% zeros(N, N);
	x       = zeros(N, N);   y = x;              bx  = x;
    by      = x;             W = x;              Wx  = ones(size(F)); Wy = Wx;
 	
	% Decompose the matrix R and G
	[Ur, Lr]= eig(transpose(R)*R);         Ur  = Ur';
	[Ug, Lg]= eig(G*transpose(G));         Ug  = Ug';
    
	% Decompose differential opertator
	%[dx dxt dy dyt] = Dxyt2(F);   
    [dx dy] = gradCal3(F, 3);
    dxt     = dx';      dyt = dy';
    
	[Ux, Lx] = eig(dxt*dx);                 Ux = Ux';
	[Uy, Ly] = eig(dy*dyt);                 Uy = Uy';
    
	% The denomitor --> recheck
	denoF = nuy + mu*diag(Lr,0)*diag(Lg,0)';
	denoV = nuy + lambda* ( diag(Lx,0)*ones(1,N) + (diag(Ly,0)*ones(1,N))' );
    denoV1 = nuy + (lambda+opts.gamma)* ( diag(Lx,0)*ones(1,N) + (diag(Ly,0)*ones(1,N))' );
    
    %% Main loop
%     eachPSNR = zeros(nOuter);       eachErr = eachPSNR;	 
%     eachSSIM = eachPSNR;
	errorKPlus = 10;                        outer = 0; 
    
    %% update parameters
    par.R = R;  par.G = G;  
    par.dx=dx;  par.dy=dy;  
    par.Ur=Ur;  par.Lr=Lr;  par.Ug=Ug;  par.Lg=Lg;
    par.Ux=Ux;  par.Lx=Lx;  par.Uy=Ly;  par.Ly=Ly;
    par.init = 0;
    par.V = V;  par.x = x;  par.y = y;  par.F = F;
    par.bx=bx;  par.by=by;  par.W = W;
    par.denoF = denoF;      par.denoV = denoV;      par.denoV1=denoV1;
    nLoop = nInner;
    
    %% main loop
    while ((errorKPlus> tol) && outer ~= nOuter)
        outer = outer +1;       F1 = F;        
        
        % initial image
        if(outer < 5)            
            [res]       = AWTV(Y, lambda, nuy, mu, nLoop, Wx, Wy, par); % 20 loop, no weight
        else
            % Weight calculation
            if(opts.weight_mode ~= 1 && outer > 150)
                 [Wx, Wy] = Weight_Assign(opts.weight_mode, par.F, par.x, par.y, opts.tau1, opts.hWeight);
            end;
            
            if(mode == 2) % spatial regularization
                nl_spatial      = post_filter(1, par.V, opts.hReg, zeros(size(Wx)));                
%                 par.x           = dx*nl_spatial;
%                 par.y           = nl_spatial*dy;
                [res]           = AWTVNL1(Y, lambda, nuy, mu, nLoop, Wx, Wy, nl_spatial, par);            
            else % gradient regularization
                if(mode == 3)
                    nl_spatial  = post_filter(1, par.V, opts.hReg, zeros(size(Wx)));
%                     nl_spatial = NLMF(par.F/(2*max(abs(par.F))), Options);  nl_spatial = nl_spatial * (2*max(abs(par.F)));
                    nl_grX      = dx * nl_spatial;
                    nl_grY      = nl_spatial * dy;
%                     par.x       = dx*nl_spatial;
%                     par.y       = nl_spatial*dy;
                elseif(mode ==4)
                    nl_grX      = post_filter(1, dx*par.V, opts.hReg, zeros(size(Wx)));
                    nl_grY      = post_filter(1, par.V*dy, opts.hReg, zeros(size(Wx)));
                    if(opts.weight_mode ~= 1 && outer > 10)
                         [Wx, Wy] = Weight_Assign_filterOut(opts.weight_mode, par.F, par.x, par.y, opts.tau1, opts.hWeight);
                    end;
                end;
                [res] = AWTVNL2(Y, lambda, nuy, mu, nLoop, Wx, Wy, nl_grX, nl_grY,par);                
            end;
        end;
        
        F = res.F;        
        par.V = res.V;  par.x = res.x;  par.y = res.y;
        par.bx=res.bx;  par.by=res.by;  par.W = res.W;
        par.F = res.F;
        
		Y = Y+Y0-R*F*G; 
        
        errorKPlus = sum((F(:) - F1(:)).^2)/sum(F(:).^2); 
        if(opts.isShowPSNR == 1)
            eachErr(outer) = errorKPlus ;
            eachPSNR(outer) =  psnr(F, image1);        
            eachSSIM(outer)= cal_ssim(F, image1, 0, 0 );           
        
            display(['outer' num2str(outer) ', inner' num2str(nInner) ', PSNR:' num2str(eachPSNR(outer)) ...
                 ', SSIM:' num2str(eachSSIM(outer)) ', errK+:' num2str(eachErr(outer))]);
        end; 
%         figure(1);
%         imshow(F, []);
    end
    %% results
    result.eachErr = errorKPlus;
    if(opts.isShowPSNR == 1)
        result.eachPSNR = eachPSNR;
        result.eachErr = eachErr;
        result.eachSSIM = eachSSIM;
    end
%     t = toc;
    
    function y = psnr(im1,im2);
        [m,n] = size(im1);
        x1 = double(im1(:));
        x2 = double(im2(:));
        mse = norm(x1-x2);
        mse = (mse*mse)/(m*n);
        if mse >0
            y = 10*log10(255^2/mse);
        else
            disp('infinite psnr');
        end
    end;

end
