function [y_dec, bits, entp, dict] = Custormized_Huffman_Coding_Laplacian(y)
%   Function for Huffman coding frame base for residual signal with uniform
%   quantization
%   - Input :  
%       + residual measurement y, follow laplacian distrinution
%   - Output:
%       + y_rec: dequantized measurement                                   
%       + bits: number of bits, should devide to total pixel to get bpp
%       + entp: entropy of measurement, should divide to toal pixel to get
%       real entropy
%       + dict: dictionary given by huffman code
    
    max_quantized_y = max(y(:));
    min_quantized_y = min(y(:));
	
	% II. Send the signal
		% create Huffman dictionary
			y_col = y(:);			
			symbols = min_quantized_y:1:max_quantized_y;
					
			prob_dist_org = hist(y_col, symbols);       
			prob_dist_org = prob_dist_org / length(y_col);
			[prob_dist, bits_for_pdf_indication] = Huff_prob_est(prob_dist_org, symbols);            
			dict = huffmandict(symbols,prob_dist); 			
			%draw_figure_entropy(prob_dist_org, prob_dist, min_quantized_y, max_quantized_y);
			
		% entropy encode  
			hcode = huffmanenco(y_col, dict); % Encode the data.
		
		% Calculate related information
			% rate
				bits = length(hcode) + bits_for_pdf_indication;
			% entropy
				entp = Measurement_Entropy_Frame(y);
				
		% entropy decode
			y_col_dec = huffmandeco(hcode,dict); % Decode the code    
            y_dec = reshape(y_col_dec, size(y));	
		
	% IV. Check the errorkid
		error = sum(sum(y_dec - y));
end

function [prob_dist, bits_for_pdf_indication] = Huff_prob_est(prob_dist_org, symbols)
    
    % notify min/max value
		min_max_bits = 64;			
    
    % transmit laplacian distribution
		alpha_bits = 64;	% 64 bits for alpha transmitting
		bits_for_pdf_indication = alpha_bits + min_max_bits;    
		alpha = 2 * max(prob_dist_org); 
		
		k = 1;
		for  i = min(symbols):1:max(symbols)
			prob_dist(k) = (alpha / 2) * 2.71828 ^ (- alpha * abs(i));
			k = k + 1;
		end
		prob_dist = prob_dist / (sum(prob_dist));
end

function draw_figure_entropy(prob_dist_org, prob_dist, min_quantized_y, max_quantized_y)
    color  = ['r', 'g', 'b', 'c', 'm', 'y', 'b', 'w'];
    marker = ['s', 'd', 'o', 'v', '*', '>', '<', '+', 'None'];
    line   = ['-', '--', ':', '-.'];

	% I. q = 8
		% II.1. blk = 16
            figure_id = 1;
			figure(figure_id); hold on;  grid on;  set(gca, 'Fontsize', 13);    
                x = min_quantized_y:1:max_quantized_y;
                bar(x, prob_dist_org);
                plot (x, prob_dist,  'Color', color(2), 'Marker', marker(9), 'Line', line(4), 'LineWidth', 3); % axis([0 0.6 28 40]);
                
				ylabel('probability'); xlabel('residual y_r_e_s');
				title(['estimated Laplacian distribution for Lenna residual @ q = 16 & blk = 16']); 
				legend( 'histogram of residual', 'estimated Laplacian distribution',  ...
						'Location', 'NorthEast');
				hold off 		
end	

