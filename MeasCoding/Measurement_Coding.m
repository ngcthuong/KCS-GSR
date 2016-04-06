function [Yq, rate, entrp]  = Measurement_Coding(Y, qStep, quan_mode_id, img_size, source_dist_type)
%   function for quantization of measurement
%   - Input:
%       + Y: input measurement
%       + qStep: Quantization step size
%       + isQuantize: quantization mode, 0 - no Quant, 1-SQ, 2-CQ

quant_mode = {'No', 'SQ'};

switch quant_mode{quan_mode_id}
    case 'No'
        Yq          = Y;
        rate        = 0;
        entrp        = 0;

    case 'SQ'
        [Yq, rate, entrp]   = SQ_Coding(Y, qStep, img_size, img_size, source_dist_type);
        
    case 2
        
end;