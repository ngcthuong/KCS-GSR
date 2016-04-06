function [rec psnrMH] = postMH(R, G, Y, initial_image, original_image, w, test_lambda, nLoopMH, ratio, optsMH)
%   function to post processing image by MH
%   Input:
%   - R, G: KCS sensing matrices and measurement Y
%   - initial_image: initial reconstructed image
%   - original_image: the original image 
%   - MH parameters: w, test_lambda, ratio
%   - optsMH: parameters for TV to reconstruct resiual image, 
%   Output:
%   - rec: reconstructed image
%   - psnrMH and ssimMH

rec = initial_image;

psnrMH = zeros(2, nLoopMH);
for iter = 1:nLoopMH
    psnrMH(1, iter) = psnr(rec, original_image);
    P_mh            = MH_Predictions_new(Y, R, G, rec, w ,test_lambda, ratio);
    psnrMH(2, iter) = psnr(P_mh, original_image);
    
    % ------ residual reconstruction -------
    [R_mh ~]        = recSepWTV(R, G, Y - R*P_mh*G, optsMH, original_image);
    rec             = R_mh + P_mh;
end;

