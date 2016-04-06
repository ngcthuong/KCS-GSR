function [F result] = recSepTV_Org(R, G, Y, opts, image1);
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
    lambda = opts.lambda;
    nuy = opts.nuy ;
    mu = opts.mu ;
    nInner = opts.nInner;
    nOuter = opts.nOuter;
    tol = opts.tol;
    
%% Initial value    
    tic
    [M, N] = size(R);       
    % Reserve memory for the auxillary variables
    Y0 = Y;             F = zeros(N,N);     V = F;			% zeros(N, N);
	x = zeros(N, N);    y = x;              bx = x;
    by = x;             W = x;
 	
	% Decompose the matrix R and G
	[Ur, Lr] = eig(transpose(R)*R);         Ur = Ur';
	[Ug, Lg] = eig(G*transpose(G));         Ug = Ug';
    
	% Decompose differential opertator
	%[dx dxt dy dyt] = Dxyt2(F);   
    [dx dy] = gradCal3(F, 3);
    dxt = dx';      dyt = dy';
    
	[Ux, Lx] = eig(dxt*dx);                 Ux = Ux';
	[Uy, Ly] = eig(dy*dyt);                 Uy = Uy';
    
	% The denomitor --> recheck
	denoF = nuy + mu*diag(Lr,0)*diag(Lg,0)';
	denoV = nuy + lambda* ( diag(Lx,0)*ones(1,N) + (diag(Ly,0)*ones(1,N))' );
    
    %% Main loop
    eachPSNR = zeros(nOuter, nInner);       eachErr = eachPSNR;	 
    eachSSIM = eachPSNR;
	errorKPlus = 10;                        outer = 0; 
    [h w] = size(F); 
    while ((errorKPlus> tol) && outer ~= nOuter)
        outer = outer +1;
        inner = 0;  
        psnrAvg = zeros(1, nInner);         errAvg = psnrAvg;
        ssimAvg = psnrAvg;
        
        Yk = mu*(R)'*Y*(G)';	
        while((inner ~= nInner)&& errorKPlus> tol) 
            inner = inner + 1;          F1 = F;         
            
            % update F   			
            Fk = Ur*(Yk + nuy*(V + W))*Ug'; 
			Fk = Fk./denoF;			
			F = (Ur)'*Fk*Ug;           
						
			% Update V			
			Vk = (Ux*(lambda* (dxt*(x-bx) + (y - by)*dyt) + nuy*(F - W))*(Uy)');
			Vk = Vk./denoV;
			V = (Ux)'*Vk*Uy;						
			
            % update x and y			
            grVx = dx*(V) + bx;
            grVy = (V)*dy + by;
            x = max(abs(grVx) - 1/lambda, 0).*sign(grVx);
            y = max(abs(grVy) - 1/lambda, 0).*sign(grVy);
            
            bx = grVx-x;
            by = grVy-y;
			W = V + W - F;            
			
            % compare results;            
            errorKPlus = norm(F - F1,2)/norm(F,2);
            errAvg(inner) = errorKPlus;
            psnrAvg(inner) = csnr(F(2:h-1,2:w-1), image1(2:h-1,2:w-1), 0, 0 );
            ssimAvg(inner) = cal_ssim(F(2:h-1,2:w-1), image1(2:h-1,2:w-1), 0, 0 );
        end       
		Y = Y+Y0-R*F*G;   
        eachPSNR(outer,:) = psnrAvg;
        eachErr(outer,:) = errAvg;  
        eachSSIM(outer, :)= ssimAvg;             
        
        display(['outer' num2str(outer) ', inner' num2str(inner) ', PSNR:' num2str(psnrAvg(inner)) ...
                 ', SSIM:' num2str(ssimAvg(inner)) ', errK+:' num2str(errorKPlus)]);
           
    end
    result.eachPSNR = eachPSNR;
    result.eachErr = eachErr;
    result.eachSSIM = eachSSIM;
    
%     t = toc;
            
     
