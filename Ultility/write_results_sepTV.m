function write_results_sepTV(file_name, subrate, psnr_final, ssim_final, rate, entr, psnr_org, ssim_org, t_org, t_post, rec_psnr_post_full)
% Function to  printf out results for SepTV package
%   write_results(file_name, subrate, psnr_final, ssim_final, rate, entr, psnr_org, ssim_org, rec_psnr_post_full)
%   Input: 
%       - file_name: file name of txt to write data
%       - subrate: current subrate
%       - psnr_final: vector include psnr of reconstructed image with
%       no_average_time (ussually 5) after post process.
%       - ssim_final: vector ssim of reconstructed image after post proces.
%       - rate: vector rate/bpp
%       - entr: vector entropy, 
%       - psnr_org: vector psrn without post process
%       - ssim_org: vector of ssim of without post process
%       - rec_psnr_post_full: the vector include psnr of several times post
%       process at last iteration. (including filtered, and filtered + res)
%   Output: 
%       - file_name.txt: include data


fid = fopen(file_name, 'a+');

fprintf(fid, '             %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f', subrate, mean(psnr_final), mean(ssim_final), ...
                    mean(rate), mean(entr), mean(psnr_org), mean(ssim_org), mean(t_org), mean(t_post));	

fprintf(fid, '      ');

% post psnr of last iteration Fitered + residual, and Filtered
fprintf(fid, ' -- ');
for i = 1:length(rec_psnr_post_full(1, :))
    fprintf(fid, '%6.3f ', rec_psnr_post_full(1, i));
end;

fprintf(fid, ' -- ');
for i = 1:length(rec_psnr_post_full(2, :))
    fprintf(fid, '%6.3f ', rec_psnr_post_full(2, i));
end;

% psnr final 5 iteration 
fprintf(fid, ' -- ');
for i = 1:length(psnr_final)
    fprintf(fid, '%6.3f ', psnr_final(i));
end;

% psnr & ssim org 5 iteration
fprintf(fid, ' -- ');
for i = 1:length(psnr_final)
    fprintf(fid, '%6.3f ', psnr_final(i));
end;

fprintf(fid, '\n');

fclose(fid);
