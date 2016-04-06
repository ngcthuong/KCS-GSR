function [rec, rec_psnr, rec_ssim, rec_psnr_full] = PostProcessing(post_mode_id, R, G, Y, rec_input, opts, image_org, imgName)

opts_post = opts; 
opts_post.nInner = opts.nInnerPost;
opts_post.nOuter = opts.nOuterPost;
post_mode = {'No', 'BM3D', 'MH', 'NLM'};
switch post_mode{post_mode_id}
    case 'MH'
        [rec, rec_psnr_full] = postMH(R, G, Y, rec_input, image_org, opts.w, opts.test_lambda, opts.nLoopPost, opts.ratio,  opts_post);
    case 'BM3D'
        [rec, rec_psnr_full] = postBM3D(R, G, Y, rec_input, image_org, opts.sigma, opts.profile, opts.nLoopPost, opts_post);
    case 'NLM'
        
    case 'No'
        rec                  = rec_input;
        rec_psnr_full        = zeros(2, opts.nLoopPost);        
end
if(opts.isShowPSNR == 1)
    rec_psnr = psnr(rec, image_org);
    rec_ssim = cal_ssim(rec, image_org, 0, 0);          
    display(['    ', imgName, ', Post Mode: ', post_mode{opts.post_mode}, ', PSNR: ', num2str(rec_psnr) 'dB, SSIM:' num2str(rec_ssim)]);
else
    rec_psnr = 0;
    rec_ssim = 0;          
end