function [res ] = SepTV_Recovery(rec_mode_id, Y, Wx, Wy, nl_spatial, nl_grx, nl_gry, opts, image1)
% function for general reconstruction of Cartoon and texture image with
% inherit the previous output
%   [res ] = SepTV_Recovery(Y, Wx, Wy, nl_spatial, nl_grx, nl_gry, opts)
%   - Input:
%       + Y: KCS measurement
%       + Wx: input weight for gradient image x
%       + Wy: input weight for gradient image y
%       + nl_spatial: filtered the spatial iamge, for NLR1, and NLR2
%       + nl_grx: filtered gradient image x
%       + nl_gry: filtered gradient image y
%   - Output:
%       + rec: the output struct include all temporal parameters such as F,
%       V, x, y, bx, by, W, 

rec_mode        = {'ATV', 'TVNL1', 'TVNL2', 'TVNL3'};

switch rec_mode{rec_mode_id}
    case 'ATV' % ATV
        [res] = AWTV(Y, opts.lambda, opts.nuy, opts.mu, opts.nLoop, Wx, Wy, opts);
        
    case 'TVNL1'
        [res] = AWTVNL1(Y, opts.lambda, opts.nuy, opts.mu, opts.nLoop, Wx, Wy, nl_spatial, opts);
        
    case 'TVNL2'
        [res] = AWTVNL2(Y, opts.lambda, opts.nuy, opts.mu, opts.nLoop, Wx, Wy, nl_grx, nl_gry, opts);
        
    case 'TVNL3' %	--- should disable the filter Cartoon   
        % 2-NLM, 3-BM3D
        nl_grx      = post_filter(2, nl_grx, opts.sigma, zeros(size(Wx)));
        nl_gry      = post_filter(2, nl_gry, opts.sigma, zeros(size(Wy)));
        [res]       = AWTVNL2(Y, opts.lambda, opts.nuy, opts.mu, opts.nLoop, Wx, Wy, nl_grx, nl_gry, opts);
        
end;
%display(['            PSNR: ' num2str(psnr(res.F, image1))]);