function [] = save_image_result(rec, folder_name, img_name, subrate, rec_psnr, rec_ssim, rec_rate)

patch = ['ResultImage\' folder_name];
if ~exist(patch, 'dir');
    mkdir(patch);
end;

final_name = [img_name, '_Sub', num2str(subrate), '_PSNR', num2str(rec_psnr), '_SSIM', num2str(rec_ssim) ];
if(rec_rate ~= 0)
    final_name = [final_name '_rate' num2str(rec_rate)];
end
final_name = [final_name '.tif'];
imwrite(uint8(rec), strcat([patch '\', final_name]));