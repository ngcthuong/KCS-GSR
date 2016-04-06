function [F, result] = DecWTVNLR(R, G, Y, opts, image1)
% Function for Decomposition based TV recovery
%   [F, result] = DecWTVNLR(R, G, Y, opts, image1)
%   Step:
%       1. Initial image is reconstructed by ATV, 
%       2. Reconstruct cartoon by opts.car.rec_mode
%       3. filtered the output reconstruction image by NLM/BM3D or similar
%       4. reconstruct the residual image by opts.tex.rec_mode
%       5. come back to step 2 and do it iteratively
%   Note: all the results from the previous iteration is used as initial of
%   the next iteration

%% test function
    % NLM parameter setup.
%     load testLenaS01
    Options.kernelratio = 6 ;%: %Radius of local Patch (default 3)
    Options.windowratio = 3;%: Radius of neighbourhood search window (default 3)
    Options.nThreads    = 4;
    Options.blocksize   = max(size(image1));
    
%     Options.enablepca   = 1;
    opts.car.Options    = Options;
    opts.tex.Options    = Options;
%% Initial value    
    tic
    [~, N] = size(R);       
    % Reserve memory for the auxillary variables
        Y0      = Y;                F       = zeros(N,N);     V       = F;			% zeros(N, N);
        x       = zeros(N, N);      y       = x;              bx      = x;
        by      = x;                W       = x;
        Wx_car  = ones(N, N);       Wy_car  = Wx_car;
        Wx_tex  = Wx_car;           Wy_tex  = Wy_car;
    
	% Decompose the matrix R and G
        [Ur, Lr] = eig(transpose(R)*R);                 Ur      = Ur';
        [Ug, Lg] = eig(G*transpose(G));                 Ug      = Ug';
    
	% Decompose differential opertator
        [dx, dy] = gradCal3(F, 3);  dxt = dx';          dyt     = dy';

        [Ux, Lx] = eig(dxt*dx);     Ux  = Ux';
        [Uy, Ly] = eig(dy*dyt);     Uy  = Uy';
    
	% The denomitor
	% cartoon
        opts.car.denoF  = opts.car.nuy + opts.car.mu*diag(Lr, 0)*diag(Lg, 0)';
        opts.car.denoV  = opts.car.nuy + opts.car.lambda*(diag( Lx, 0) * ones(1, N) + (diag( Ly, 0) * ones(1, N))');
        opts.car.denoV1 = opts.car.nuy + (opts.car.lambda + opts.car.gamma)*(diag( Lx, 0) * ones(1, N) + (diag( Ly, 0) * ones(1, N))');
    % texture
        opts.tex.denoF  = opts.tex.nuy + opts.tex.mu*diag(Lr, 0)*diag(Lg, 0)';
        opts.tex.denoV  = opts.tex.nuy + opts.tex.lambda*(diag( Lx, 0) * ones(1, N) + (diag( Ly, 0) * ones(1, N))');
        opts.tex.denoV1 = opts.tex.nuy + (opts.tex.lambda + opts.tex.gamma)*(diag( Lx, 0) * ones(1, N) + (diag( Ly, 0) * ones(1, N))');
 
    
    %% Updating parameters --> check the  
	% for Cartoon
    opts.car.R = R;  opts.car.G = G;  
    opts.car.dx=dx;  opts.car.dy=dy;  
    opts.car.Ur=Ur;  opts.car.Lr=Lr;  opts.car.Ug=Ug;  opts.car.Lg=Lg;
    opts.car.Ux=Ux;  opts.car.Lx=Lx;  opts.car.Uy=Ly;  opts.car.Ly=Ly;  
    
    opts.car.init = 0;
    opts.car.V = V;  opts.car.x = x;  opts.car.y = y;  opts.car.F = F;
    opts.car.bx=bx;  opts.car.by=by;  opts.car.W = W;
	
	% -- Texture
    opts.tex.R = R;  opts.tex.G = G;  opts.tex.dx=dx;  opts.tex.dy=dy;  
    opts.tex.Ur=Ur;  opts.tex.Lr=Lr;  opts.tex.Ug=Ug;  opts.tex.Lg=Lg;
    opts.tex.Ux=Ux;  opts.tex.Lx=Lx;  opts.tex.Uy=Ly;  opts.tex.Ly=Ly;  
    
    opts.tex.init = 0;
    opts.tex.V = V;  opts.tex.x = x;  opts.tex.y = y;  opts.tex.F = F;
    opts.tex.bx=bx;  opts.tex.by=by;  opts.tex.W = W;
	
    opts.tex2 = opts.tex; 
	
    % init measurement
    Y_tex_hat   = 0.*Y0;
    
    %% main loop
    
    eachPSNR    = zeros(opts.nOuter);        eachErr = eachPSNR;	     eachSSIM = eachPSNR;
	errorKPlus  = 10;                   outer   = 0; 
    
    while ((errorKPlus> opts.tol) && outer ~= opts.nOuter)
        outer = outer +1;       F1 = F;        
        
        %% reconstruction mode
		if(opts.dtv_mode == 0) % residual mode
			Y_car = Y;
		elseif (opts.dtv_mode == 1) % Cartoon texture mode
			Y_car = Y - Y_tex_hat;
		end;
		
        %% cartoon image
        display(['Outer: ' num2str(outer)]);
        if(outer == 1) %initial iteration
            [res_car] = AWTV(Y_car, opts.car.lambda, opts.car.nuy, opts.car.mu, opts.car.nLoop_init, Wx_car, Wy_car, opts.car);            
        else			
			% weight calculation
            % Take from filtered final image, 
            tmpX    = dx*F;     tmpY = F*dy; 
            %opts.car.x = tmpX;  opts.car.y = tmpY;
            [Wx_car, Wy_car] = Weight_Assign_filterOut(opts.car.weight_mode, F, tmpX, tmpY, opts.car.tau, opts.car.sigma);
            
			% recovery the Cartoon (main part)
            % res_car.F is the F from previous iteration --->  
            
           % [res_car] = SepTV_Recovery(opts.car.rec_mode, Y_car, Wx_car, Wy_car, res_car.F, ...
            %                           dx*res_car.F, res_car.F*dy, opts.car, image1);
			[res_car] = SepTV_Recovery(opts.car.rec_mode, Y_car, Wx_car, Wy_car, F, ...
                                      dx*F, F*dy, opts.car, image1);
        end;
        
		tmp_car       = res_car.F;
        
        % filter the Cartoon image -- Optional
        res_car.F   = post_filter(opts.car.filter_mode, res_car.F, 35, image1);                

        opts.car.V  = res_car.V;    opts.car.x      = res_car.x;  
        opts.car.y  = res_car.y;   
        
        % Cartoon measurement
        Y_car_hat   = R*res_car.F*G;     
         
        %% Texture/residual Image
        % Input texture measurement
        Y_tex       = Y - Y_car_hat;
        
%         display(['Texture: ']);
        if(outer == 1)
            opts.tex.nLoop_init    = 10;
            [res_tex]        = AWTV(Y_tex, opts.tex.lambda, opts.tex.nuy, opts.tex.mu, opts.tex.nLoop_init, Wx_tex, Wy_tex, opts.tex);
            res_tex2 = res_tex; 
        else
			% weight matrix for texture part
%             [Wx_tex, Wy_tex] = Weight_Assign(opts.tex.weight_mode, res_tex.F,  dx*res_tex.F, res_tex.F*dy, opts.tex.tau, opts.tex.sigma);
            
			% recovery the Cartoon (main part)
            [res_tex]        = SepTV_Recovery(opts.tex.rec_mode, Y_tex, Wx_tex, Wy_tex, res_tex.F, ...
                                       dx*res_tex.F, res_tex.F*dy, opts.tex, image1);
            
        end;
		
		%% filter the Textures image --Optional
        tmp_tex     = res_tex.F;
		% for Decomposition only, not as good as filter combined image
% 		res_tex.F   = post_filter(2, res_tex.F, 0.0025, image1);        
        
        opts.tex.V  = res_tex.V;    opts.tex.x      = res_tex.x;  
        opts.tex.y  = res_tex.y;   
        
        %% Texture measurement
        Y_tex_hat   = R*res_tex.F*G;    
        
        %% Combine image
        F_tmp       = res_car.F + res_tex.F;   
        F1           = post_filter(opts.tex.filter_mode, F_tmp, 15, image1);   
       
        %% the third layer
        Y_tex       =  Y - R*F1*G;
        [res_tex2]  = SepTV_Recovery(opts.tex2.rec_mode, Y_tex, Wx_tex, Wy_tex, res_tex2.F, ...
                                       dx*res_tex2.F, res_tex2.F*dy, opts.tex2, image1);
        F_tmp       = F1 + res_tex2.F;   
        F2          = post_filter(opts.tex.filter_mode, F_tmp, 7, image1);   
                       
        
        F           = F2; 
        %% update measurement
		Y           = Y+Y0-R*F*G; 
        
        %% error 
        errorKPlus = sqrt(sum((F(:) - F1(:)).^2)/sum((F(:).^2))); 
        eachErr(outer) = errorKPlus ;
        eachPSNR(outer) =  psnr(F, image1 );   
        
        %% show psnr 
            % Cartoon only  
            psnr1(outer) = psnr(res_car.F, image1); 
            psnr2(outer) = psnr(F1, image1); 
            psnr3(outer) = psnr(F2, image1); 
            res_car.F    = F2;
            display(['     Cartoon,: ' num2str(psnr1(outer)) ', F1:' num2str(psnr2(outer)) ',F2:' num2str(psnr3(outer))]);
%             psrn_car_woFilter(outer) = psnr(tmp_car, image1);
%             psrn_car_wFilter(outer)  = psnr(res_car.F, image1);
%             display(['     Cartoon, woFilter: ' num2str(psrn_car_woFilter(outer)) ', wFilter:' num2str(psrn_car_wFilter(outer))]);
%             
%             % Combine 
%             psrn_com_woFilter(outer) = psnr(F_tmp, image1);
%             psrn_com_wFilter(outer)  = psnr(F, image1);
%             psrn_com_wFilter_mix(outer)  = psnr(res_car.F + res_tex.F, image1);
%             display(['     Comb, woFilter: ' num2str(psrn_com_woFilter(outer)) ', wFilter:' num2str(psrn_com_wFilter(outer)) ...
%                      ', wFilter mix:' num2str(psrn_com_wFilter_mix(outer))]);
             
%             figure(1); subplot(1,2,1); imshow(tmp_car, []); title(['Car iter' num2str(outer) ', PSNR:' num2str(psrn_car_woFilter(outer))]);
%             figure(1); subplot(1,2,2); imshow(res_car.F, []); title(['Car filtered iter ' num2str(outer) ', PSNR:' num2str(psrn_car_wFilter(outer))]);
%             figure(2); subplot(1,2,1); imshow(tmp_tex, []); title(['Tex iter' num2str(outer) ]);
%             figure(2); subplot(1,2,2); imshow(res_tex.F, []); title(['Tex filtered iter' num2str(outer) ]);
%             figure(3); subplot(1,2,1); imshow(F_tmp, []); title(['Comb iter' num2str(outer) ', PSNR:' num2str(psrn_com_woFilter(outer))]);            
%             figure(3); subplot(1,2,2); imshow(F, []); title(['Comb Filtered iter' num2str(outer) ', PSNR:' num2str(psrn_com_wFilter(outer))]);            
%             pause
        display(['outer' num2str(outer) ', PSNR:' num2str(eachPSNR(outer)) ...
                 ', errK+:' num2str(eachErr(outer))]);

    end 
    result.eachPSNR = eachPSNR;
    result.eachErr = eachErr;
%     result.eachSSIM = eachSSIM;
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

