function H1 = MH_Zero_Padding_new(cur_block, i, j, rows, cols, ratio)
% H = MH_Zero_Padding_new(cur_block, i, j, rows, cols, ratio)
%	Function to padding zero into the current block 
%   Input: 
%       + cur_block: the current block, size 64x64
%       + i, j: location index of current block, i(1:4), j(1:4)
%       + rows, cols: size of the current image , 256x256
%       + ratio: number of block, r = 4, total blocks 4x4 = 16

% for ratio = 4
blk_size = rows/ratio;
H = zeros(rows, cols);
H( (i-1)*blk_size + 1:   i*blk_size, (j - 1)*blk_size +1: j * blk_size) = cur_block;
H1 = H(:);
