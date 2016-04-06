function predicted_image = MH_Predictions_new(Y, R, G, reference, w ,test_lambda, ratio)
% predicted_image = MH_Predictions_new(Y, R, G, reference, w ,test_lambda, ratio)
%	Function to generate predicted image using Multi-hypothesis
%	- Input: 
%		+ Y: received measurement, size m x m
%		+ R, G: Kronecker Sensing Matrix, size m x n,
%		+ reference: reference image for MH searching, size n x n
%		+ w: extended search range; window size (current_size + w), example (10-20);
%		+ test_lambda: labmda value for Tikhonov regularization; example 1.0 
%		+ ratio: the ratio for small block size, example: 4, for split image into 4x4 = 16 block
%	- Output:
%		+ predicted_image:

% load testMHLenaS02

[rows cols] = size(reference);
[yrows ycols] = size(Y);
blk_size = rows / ratio;

predicted_image = zeros(rows, cols);

% extend image by mirroring
% ext_ref = symextend(reference, w);
ext_ref = padarray(reference,[w w], 'symmetric', 'both');

% Sensing function
P = @(x) R * x * G;

% extend window
% number of reference
noRef = (w*2 + 1)^2;
H1 = zeros(rows * cols, ratio^2 * noRef);
A  = zeros(yrows*ycols, ratio^2 * noRef);

index = 0;
for i = 1: ratio
	for j = 1: ratio
		x = (i-1)*blk_size + w; y = (j-1)*blk_size + w;
		
		ref_blk = ext_ref( x - w +1: x + blk_size + w,  y - w +1: y + blk_size + w);
        tmpBlk = im2col(ref_blk, [blk_size, blk_size], 'sliding');
        
        for k = 1: size(tmpBlk, 2)
            index = index + 1;
            curBlk = reshape(tmpBlk(:, k), [blk_size blk_size]);
            H1(:, index) = MH_Zero_Padding_new(curBlk, i, j, rows, cols, ratio);
            tmpA = R*reshape(H1(:, index), rows, cols)*G;
            A(:, index) = tmpA(:);
        end;
	end;
end;
cur_proj = Y(:);
% A = P(H); % projected hypotheses matrix
norms = A - repmat(cur_proj,[1 size(A,2)]);
norms = sum(norms.^2);

lambda = adaptive_lambda(cur_proj, A);

G = diag(norms);     % diagonal Tikhonov matrix
weights = (A'*A + lambda.*G)\(A'*cur_proj);
predicted_image = reshape(H1 * weights(:), [rows cols]);
