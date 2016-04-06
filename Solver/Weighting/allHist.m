function [Grad Energy TV hisGrad] = allHist(gradX, gradY, nbrHis, step)
nbrHis = floor(nbrHis./step);
hisGrad = zeros(2, nbrHis);
    
    %hisGrad(1, 1) = nnz(abs(gradX) < step);
    %hisGrad(2, 1) = nnz(abs(gradY) < step);
    tmpGradX = 0*gradX;
    tmpGradY = 0*gradY;
    
    for i = 1 : nbrHis
        hisGrad(1, i) = nnz(abs(gradX) < step*i);% - hisGrad(1, i-1);
        hisGrad(2, i) = nnz(abs(gradY) < step*i); %- hisGrad(2, i-1);
        
        tmp1 = (abs(gradX) < i*step) - tmpGradX;
        tmpGradX = tmpGradX + tmp1;
        
        tmp2 = (abs(gradY) < i*step) - tmpGradY;
        tmpGradY = tmpGradY + tmp2;
        
        % Energy:
        tmpEnergy = (tmp1.*gradX);
        Energy(1,i) = sum(tmpEnergy(:).^2);
        TV(1,i) = sum(abs(tmpEnergy(:)));
        
        tmpEnergy = (tmp2.*gradY);
        Energy(2,i) = sum(tmpEnergy(:).^2);
        TV(2,i) = sum(abs(tmpEnergy(:)));
    end
    Grad = hisGrad;
    Grad(:, 2:nbrHis) = hisGrad(:, 2:nbrHis) - hisGrad(:, 1:(nbrHis-1));
%     Grad(1,:) = Grad(1,:)/(hisGrad(1, nbrHis));
%     Grad(2,:) = Grad(1,:)/(hisGrad(2, nbrHis));
    
    Grad(1,:) = Grad(1,:)/ max(Grad(1,:));
    Grad(2,:) = Grad(2,:)/ max(Grad(2,:));
    TV(1,:) = TV(1,:)/ max(TV(1,:));
    TV(2,:) = TV(2,:)/ max(TV(2,:));