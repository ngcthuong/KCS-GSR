function [rec, rec_psnr, rec_ssim, results]         = recSep(rec_mode, R, G, Y, opts, image, imgName)

switch rec_mode
    case 1  % Org SepTV
        [rec, results]         = recSepWTV(R, G, Y, opts, image);                % Org TV
    case 2        
        [rec, results]         = recSepWTVNL(R, G, Y, opts, image, rec_mode);    % TVNLR1
    case 3
        [rec, results]         = recSepWTVNL(R, G, Y, opts, image, rec_mode);    % TVNLR2
    case 4
        [rec, results]         = recSepWTVNL(R, G, Y, opts, image, rec_mode);    % TVNLR3
    case 5
        [rec, results]         = DecWTVNLR(R, G, Y, opts, image);                % DTV
    case 6
        [opts.initial, ~]      = DecWTVNLR(R, G, Y, opts, image);
        [rec, results]         = KCS_GSR_V02(R, G, Y, opts, image);                % DTV
end;
rec_psnr = psnr(rec, image);
rec_ssim = cal_ssim(rec, image, 0, 0);
display(['    ', imgName, ', Init Rec. SepTV, PSNR:', num2str(rec_psnr), 'db, SSIM: ' num2str(rec_ssim)]);