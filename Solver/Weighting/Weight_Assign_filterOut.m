function [Wx, Wy, rec_filtered] = Weight_Assign_filterOut(weight_mode_id, rec, rec_dx, rec_dy, tau1, sigma)
% Function to calculate the weights for total variaiotni
%   [Wx, Wy] = Weight_Assign(weight_mode, rec_spatiall, rec_dx, rec_dy, tau1, tau2)
%   - Input:
%       + weight_mode: NO, ENOW, HENOW, EDGE-CS, CANDES
%       + rec_dx     : input gradient image dx
%       + rec_dy     : input gradient image dy
%       + tau1       : threshold 1; which is used in ENOW, HENOW
%       + tau2       : threshold 2; which is used in EdgeCS
weight_mode     = {'NO', 'ENOW', 'HENOW', 'edgeCS', 'Candes' };
[dx, dy] = gradCal3(rec, 3);
switch weight_mode{weight_mode_id}
    case 'NO'
        Wx      = ones(size(rec_dx));
        Wy      = Wx; 
        
    case 'ENOW'
%         rec_dx = post_filter(2, rec_dx, sigma, zeros(size(rec_dx)));
%         rec_dy = post_filter(2, rec_dy, sigma, zeros(size(rec_dx)));
%         rec     = post_filter(2, rec, sigma, zeros(size(rec_dx)));
%         rec_dx  = dx*rec;       rec_dy = rec*dy;
        Wx      = abs(rec_dx) < tau1;
		Wy      = abs(rec_dy) < tau1;
        
    case 'HENOW'        
%         rec_dx = post_filter(2, rec_dx, 0.15, zeros(size(rec_dx)));
%         rec_dy = post_filter(2, rec_dy, 0.15, zeros(size(rec_dx)));
%         rec     = post_filter(2, rec, sigma, zeros(size(rec_dx)));
%         rec_dx  = dx*rec;       rec_dy = rec*dy;
        [Grad, ~, ~, ~]   = allHist(rec_dx, rec_dy, 255, tau1); 
		[Wx, Wy]          = calWeight(rec_dx, rec_dy, Grad, tau1); 
        
    case 'edgeCS'
        
    case 'Candes'
end;
rec_filtered  = rec;