function [Gx Gy] = calWeight(grX, grY, hist, step)
[m,n] = size(grX);
Gx = zeros(m, n);    Gy = Gx;
% gLength = floor(maxGrad/step);
% gradStart = 0;
% gradEnd = step;

for i = 1: length(hist)
    gradStart = (i-1)*step;
    gradEnd = i * step;
    tmpX = hist(1, i).* ((abs(grX) < gradEnd) - (abs(grX) < gradStart));
    tmpY = hist(2, i).* ((abs(grY) < gradEnd) - (abs(grY) < gradStart));
    
    Gx = Gx + tmpX;
    Gy = Gy + tmpY;    
end
        
end