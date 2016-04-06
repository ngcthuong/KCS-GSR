% 
% function [y_q rate] = SQ_Coding(y, quantizer_bitdepth, num_rows, num_cols)
% 
%   This function performs SQ coding on the block measurements. For 
%   simplicity, this function assumes that coding and decoding occurs
%   at the same place. Therefore, the final output is recovered measurement.
%   Unlike BCS-SPL and MH-BCS-SPL, MS-BCS-SPL has different data structure. 
%   Therefore, the procedure is divided by two and handles the data 
%   accordingly.
% 
%   Inputs
%   y: block-by-block sampled measurements
%   quantizer_bitdepth: quantizer bit depth
%   num_rows: number of rows of the original image
%   num_cols: number of columns of the original image
%
%   Outputs
%   y_q: quantized measurements
%   rate: rate information reprensted as entropy

%
% DPCM for BCS-SPL: Differential Pulse-Code Modulation for 
% Block Compressed Sensing - Smooth Projected Landweber
% Copyright (C) 2012  James E. Fowler
% http://www.ece.mstate.edu/~fowler
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

%
%   Originally written by SungKwang Mun, Mississippi State University
%   Edited by Thuong Nguyen 2014 Apr. 02, Sungkyunkwan Univ.

function [y_quant, rate, entr] = SQ_Coding(y, quantizer_bitdepth, num_rows, num_cols, source_dist_type)
    y_max = max(y(:)); y_min = min(y(:));
    % q: stepsize
    q = (y_max - y_min)/2^quantizer_bitdepth;

    % simple scalar quantization
    yq = round(y/q);

	% Huffman coding
%     [y_dec, bits, entp, dict] = Custormized_Huffman_Coding_Laplacian(yq);
    [y_dec, bits, entp, dict] = Custormized_Huffman_Coding(yq, source_dist_type);
	% bit rate and entropy calculation
    total_pixels = num_rows*num_cols;
    rate = bits/total_pixels;
	entr = entp/total_pixels;
	
    y_quant = y_dec*q;
    
end