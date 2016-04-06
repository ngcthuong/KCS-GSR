function [rec] = post_filter(filter_mode_id, rec_input, sigma, image1)
% function perform diferent type of filter
%   [rec] = post_filter(filter_mode, opts)
%   - Input:
%       + filter_mode: the filter mode, BM3D, NLM or Wiener, Gaussian
filter_mode     = {'No', 'NLM', 'BM3D', 'WNNM', 'MH' , 'Wiener', 'Gauss'};
switch filter_mode{filter_mode_id}
    case 'No'
        rec = rec_input;
        
    case 'NLM'
        Options.kernelratio = 6 ;%: %Radius of local Patch (default 3)
        Options.windowratio = 3;%: Radius of neighbourhood search window (default 3)
        Options.filterstrength = sigma;
        rec         = NLMF(rec_input/255, Options);
        rec         = rec*255;
        
    case 'BM3D' % only valid for spatial domain
        [~, rec] = BM3D(image1, rec_input, sigma, 'np', 0);
        rec = rec*255;
		
    case 'WNNM' % only for spatial domain
		 Par   = ParSet(sigma);   
		 rec   = WNNM_DeNoising( rec_input, image1, Par );
	
	case 'MH' % for MH, canot use at this time
		rec = MH_Predictions_new(opts.Y, opts.R, opts.G, rec_input, opts.w , opts.test_lambda, opts.ratio);
		
    case 'Wiener'
        rec = wiener2(rec_input, [3, 3]);
        
    case 'Gauss'
        G = fspecial('gaussian',[3 3], sigma);
        rec = imfilter(rec_input, G, 'same');
end

% display(['            PostFilter' filter_mode{filter_mode_id} ', PSNR: ' num2str(psnr(rec, image1))]);
% 