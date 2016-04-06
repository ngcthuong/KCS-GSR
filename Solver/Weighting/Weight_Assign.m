function [Wx, Wy, rec_filtered] = Weight_Assign(weight_mode_id, rec, rec_dx, rec_dy, tau1, sigma)
% Function to calculate the weights for total variaiotni
%   [Wx, Wy] = Weight_Assign(weight_mode, rec_spatiall, rec_dx, rec_dy, tau1, tau2)
%   - Input:
%       + weight_mode: NO, ENOW, HENOW, EDGE-CS, CANDES
%       + rec_dx     : input gradient image dx
%       + rec_dy     : input gradient image dy
%       + tau1       : threshold 1; which is used in ENOW, HENOW
%       + tau2       : threshold 2; which is used in EdgeCS
weight_mode     = {'No', 'ENOW', 'HENOW', 'EdgeCS', 'Candes', 'PC', 'PC_ENOW', 'PC_HENOW'};
[dx, dy] = gradCal3(rec, 3);
switch weight_mode{weight_mode_id}
    case 'No'
        Wx      = ones(size(rec_dx));
        Wy      = Wx; 
        
    case 'ENOW'
%         rec_dx = post_filter(2, rec_dx, sigma, zeros(size(rec_dx)));
%         rec_dy = post_filter(2, rec_dy, sigma, zeros(size(rec_dx)));
        rec     = post_filter(2, rec, sigma, zeros(size(rec_dx)));
        rec_dx  = dx*rec;       rec_dy = rec*dy;
        Wx      = abs(rec_dx) < tau1;
		Wy      = abs(rec_dy) < tau1;
        
    case 'HENOW'        
%         rec_dx = post_filter(1, rec_dx, sigma, zeros(size(rec_dx)));
%         rec_dy = post_filter(1, rec_dy, sigma, zeros(size(rec_dx)));
        rec     = post_filter(2, rec, sigma, zeros(size(rec_dx)));
        rec_dx  = dx*rec;       rec_dy = rec*dy;
        [Grad, ~, ~, ~]   = allHist(rec_dx, rec_dy, 255, tau1); 
		[Wx, Wy]          = calWeight(rec_dx, rec_dy, Grad, tau1); 
        
    case 'edgeCS'
        
    case 'PC'
        Wx = 1 - phasecong2(rec_dx);
        Wy = 1 - phasecong2(rec_dy);
        
    case 'PC_ENOW'
        Gx1 = abs(rec_dx) < tau1;
        Gy1 = abs(rec_dy) < tau1;
        Gx2 = 1 - phasecong2(rec_dx);
        Gy2 = 1 - phasecong2(rec_dy);
        Wx  = Gx1 .* Gx2; Wy = Gy1.*Gy2;
        
    case 'PC_HENOW'
        [Grad , ~ , ~, ~ ] = allHist(rec_dx, rec_dy, 255, tau1); 
        [Gx1, Gy1] = calWeight(rec_dx, rec_dy, Grad, tau1);
        Gx2 = 1 - phasecong2(rec_dx);
        Gy2 = 1 - phasecong2(rec_dy);
        Wx  = Gx1 .* Gx2; Wy = Gy1.*Gy2;        
        
end;
rec_filtered  = rec;