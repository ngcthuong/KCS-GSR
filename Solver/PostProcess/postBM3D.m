function [rec psnrBM3D ] = postBM3D(R, G, Y, initial_image, original_image, sigma, profile, nLoopBM3D, optsBM3D)
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

psnrBM3D = zeros(2, nLoopBM3D);
for iter = 1:nLoopBM3D
    psnrBM3D(1, iter)   = psnr(rec, original_image);
    [tmpxx P_mh]        = BM3D(original_image, rec, sigma, profile, 0); P_mh = P_mh*255;
    psnrBM3D(2, iter)   = psnr(P_mh, original_image);
    
    % ------ residual reconstruction -------
    [R_mh ~]            = recSepWTV(R, G, Y - R*P_mh*G, optsBM3D, original_image);
    rec                 = R_mh + P_mh;
    
end;


