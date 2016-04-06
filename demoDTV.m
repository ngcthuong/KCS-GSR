% TV reconstruction for Kronecker Compressive Sensing
% Including original TV [1], TVNL1, TVNL2 and TVNL3
% with optional for post processing by BM3D[], Multi-hypothesis;
% [1]"S.L Shishkin, H. Wang and G.S Hagen, "Total Variation Minimization wih Separable Sensing Operator", ICISP, 2010.
% Implemented by Thuong Nguyen Canh, Feb. 2013
% Updated March 2014

% Parameter setting
% Input:
%   - N: size of input image, input image shoud be squared image
%   - sparsity test level 0.1-0.9
%   - Mode: 0 - Original TV, 1- TV+NLR1, 2-TV+NL2; 3-TV+NLR3
% For general TV
%   - opts: the general setting parameter, contains all sub-setting
%       + lamba, nuy, mu, is parameters for TV reconstruction (let as suggestion)
%       + nInner/nOuter: number of maximum outer and inner loop
%       + tol: tolerance to stop norm(F_(k+1) - F_k, 2) / norm(F_k) < tol
%       + nbrLoop: number of average times
% For regularization: need parameter for NLM filter (degree of filter)
%   - opts.hFilter: degree of filter ( 0.02 - 0.05)
%   - opts.gamma: paramter for optimization problem of TV_nonlocal
% For Post processing with MH
%   - isPostMH: with (1) or without (0) post processing with MH
%   - nLoopMH: number of time for post procesisng by MH
%   - MH parameters:
%       + w: search range extend
%       + ratio: ratio between block for MH and frame size, ie 2 or 4 or 8
%       + test_lambda: the initial lambda parameter for Tikhonov
%       regularization of MH
% For Post processing with BM3D/NLM
%   - isPostBM3D: with (1) or without (0) post processing with BM3D
%   - nLoopBM3D: number of time for post procesisng by MH

% function testSepTV_all()
function demo
close all;  clear all; % clc;
path(path,genpath(pwd));
set(0,'RecursionLimit', 1200)
%% General setting
N               = 512;
sparsity        = [0.05 0.1 0.2 0.3 0.4];
testIm          = [1 ];
% testIm          = [5 6 7 8 ];
% testIm          = [9 10 11 12 13];
rec_mode        = {'TV', 'TVNL1', 'TVNL2', 'TVNL3', 'DTV', 'GSR'} ;
weight_mode     = {'No', 'ENOW', 'HENOW', 'EdgeCS', 'Candes', 'PC', 'PC_ENOW', 'PC_HENOW'};
post_mode       = {'No', 'BM3D', 'MH', 'NLM'};
meas_dist       = {'Gauss', 'Laplacian'};
quant_mode      = {'No', 'SQ'};

rec_mode_id     = 6;
weight_mode_id  = 1;
post_mode_id    = 1;
quant_mode_id   = 1;
isShowPSNR      = 1;
nSNR            = 0;


[opts, note_str]= setup_parameter(rec_mode_id, weight_mode_id, post_mode_id, quant_mode_id, isShowPSNR);
% note_str = [ note_str '_SNR' num2str(nSNR) ];
% opts.filter_mode = 2; opts.sigma = 0.03; % NLM
opts.nbrLoop    = 5;
opts.par.ArrayNo     = 40; 
note_str = [note_str '_arrNo' num2str(opts.par.ArrayNo)];
% Quantization
quant_bit       = [0];  
% quant_bit       = [8 10];     %quant_bit_depth  = ones(size(sparsity))*quant_bit;

for mmm = 1:1:length(testIm)
    [image, img_name] = testImage(N, testIm(mmm));
    
    save_folder_text    = ['Result_text' num2str(N) '\' ];
    if ~exist(save_folder_text, 'dir');
        mkdir(save_folder_text);
    end;
    file_name_save  	= [save_folder_text 'results_' img_name '_20150323_' note_str ];
    
    if(rec_mode_id == 5)
    write_info([file_name_save  '.txt'], [img_name '_Car_lam' num2str(opts.car.lambda) '_nuy' num2str(opts.car.nuy)...
                '_mu' num2str(opts.car.mu) '_nLoopInit' num2str(opts.car.nLoop_init) '_nLoop' num2str(opts.car.nLoop)...
                '_sigma' num2str(opts.car.sigma)]);
    write_info([file_name_save  '.txt'], [ '        _Tex_lam' num2str(opts.tex.lambda) '_nuy' num2str(opts.tex.nuy)...
                '_mu' num2str(opts.tex.mu) '_nLoopInit' num2str(opts.tex.nbrLoop_init) '_nLoop' num2str(opts.tex.nLoop)...
                '_sigma' num2str(opts.tex.sigma)]);         
    else
        write_info([file_name_save  '.txt'], [img_name note_str]);
    end;
    
    for qId = 1:1:length(quant_bit)
        quant_bit_depth = quant_bit(qId);
        write_info([file_name_save  '.txt'], ['      qD:' num2str(quant_bit_depth)]);            
        Note = [note_str '_qDepth' num2str(quant_bit_depth)];        display(Note);
                
        % build the sampling matrix, R
        t_org      = zeros(opts.nbrLoop, length(sparsity));      t_post = t_org;
        psnr_final = t_org;         ssim_final  = t_org;
        psnr_org   = t_org;         ssim_org    = t_org;
        rate       = t_org;         entr        = t_org;
        rec_psnr_post = cell(1);
        
        for kk = 1:1:length(sparsity)
            
            for trial = 1:opts.nbrLoop
                NoteBeg     = ['In' num2str(opts.nInner) '_Out' num2str(opts.nOuter)];
                display([NoteBeg '_' img_name '_Sub' num2str(sparsity(kk)) '_trial' num2str(trial) '/' num2str(opts.nbrLoop )])
                
                % -------------- Sensing matrix ------------------------
                [R, G]              = KCS_SensingMtx(N, sparsity(kk), trial);
                Y                   = R*image*G;
                if(nSNR ~=0)
                noise = rand(size(Y)); %%gives a 10 by 10 matrix. 
                scale = ( std(Y(:))/std(noise(:)) ) / nSNR; 
                Y     = Y + scale * noise; 
                end;
                % -------------- Measurement quant ----------------------
                source_dist_type = 2;
                [Yq,  rate(trial, kk), entr(trial,kk)]= Measurement_Coding(Y, quant_bit_depth, opts.quant_mode, N, source_dist_type);
                Y                = Yq;
                
                display(['     Rate Total: ' num2str(rate(trial, kk)) ', Entropy: ' num2str(entr(trial, kk)) ', Loss:' ...
                                num2str( 100*(rate(trial, kk) - entr(trial, kk))/ (entr(trial, kk) + 1e-10)) '%']);
                
                % --------------- Recover the image -----------------------
                t0 = cputime;
                [rec, ~, ~, ~] = recSep(opts.rec_mode, R, G, Y, opts, image, img_name); 
                t_org(trial, kk) = cputime - t0;
                rec_org = rec;
                
                % Post pocessing
                t0 = cputime;
                [rec, ~, ~, rec_psnr_post{trial}{kk}] = PostProcessing(opts.post_mode, R, G, Y, rec, opts, image, img_name);                
                t_post(trial, kk)  = cputime - t0;
                t_post(trial, kk)  = t_post(trial, kk) + t_org(trial,kk);
                rec_post           = rec;
                
                % --------------- Final Results  ------------------------
                psnr_org(trial, kk)   = psnr(rec_org, image);
                ssim_org(trial, kk)   = cal_ssim(rec_org, image, 0, 0);
                psnr_final(trial, kk) = psnr(rec_post, image);
                ssim_final(trial, kk) = cal_ssim(rec_post, image, 0, 0);            
                
            end; % end trial
            write_results_sepTV([file_name_save '.txt'], sparsity(kk), psnr_final(:, kk), ...
                ssim_final(:, kk), rate(:, kk), entr(:, kk), psnr_org(:, kk), ssim_org(:, kk), ...
                t_org(:, kk), t_post(:, kk), rec_psnr_post{trial}{kk} );
            
            %% ============ save image  ==================
            if(opts.post_mode ~= 1)
                save_image_result(rec_org, [Note '_Org'], img_name, sparsity(kk), mean(psnr_org(:, kk)), mean(ssim_org(:)), mean(rate(:, kk)));
            end;
            save_image_result(rec_post, [ Note ] , img_name, sparsity(kk), mean(psnr_final(:, kk)), mean(ssim_final(:)), mean(rate(:, kk)));
        end; % end sparsity        
        % 	save all images
        patch3 = ['Results\' 'All_In' num2str(opts.nInner) Note];
        if ~exist(patch3, 'dir');
            mkdir(patch3);
        end;
        patch3 = [patch3 '\' img_name '_Sub' num2str(sparsity(kk)) '.mat'];
        save(patch3, 'psnr_final', 'ssim_final','psnr_org', 'ssim_org', 'rec_psnr_post' , 't_org', 't_post', 'opts', 'sparsity', 'rec') ;
        
    end; % end test image
end; % quant-step
display('END SIMULATION!!!');

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
end

end