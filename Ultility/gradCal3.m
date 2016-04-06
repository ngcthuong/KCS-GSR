function [dx dy] = gradCal3(X , type)
% [Dh Dv] = gradCal(X)
%   Function to calculate the gradient Dy is dx matrix; Dx 
%   The matrix cal is not correct for diag gradient case

    [M,N] = size(X);  m = M;
    X1 = eye(M, N);
    % The horizontal and vertial Gradient, zero version    
    switch type
        case 1 
%             [Dx Dy] = gradZero(X);
            [dx dy] = gradZero(X1);
        case 2
%             [Dx Dy] = gradCopy(X);
            [dx dy] = gradCopy(X1);
        otherwise 
%             [Dx Dy] = gradRound(X);
            [dx dy] = gradRound(X1);
    end;
             
    % The diagonal gradient
%     mtx1 = eye(m-1, m-1);
%     mtx1 = [mtx1; zeros(1, m-1)];
%     mtx1 = [mtx1, zeros(m, 1)];
%     
%     
%     mty1 = mtx1;
%     mty2 = mty1;
%         
%     % Calculate Dx2 = mtx2 * Dx
%     mtx2 = eye(m-1, m-1);
%     mtx2 = [mtx2; zeros(1, m-1)];
%     mtx2 = [zeros(m, 1), mtx2];
%     
%     outOpt.mtx1 = mtx1;
%     outOpt.mtx2 = mtx2;
%     outOpt.mty1 = mty1;
%     outOpt.mty2 = mty2;
    
end

    function [Dx Dy] = gradZero(X)
        [M,N] = size(X);
        Dx = -diff(X,[],1);
        Dx = [Dx;zeros(1,N)];
        
        Dy = -diff(X,[],2);
        Dy = [Dy zeros(M,1)];
    end
    
    function [Dx Dy] = gradCopy(X)
        [M,N] = size(X);
        Dx = -diff(X,[],1);
        Dx = [Dx;Dx(N-1, :)];
        
        Dy = -diff(X,[],2);
        Dy = [Dy Dy(:, N-1)];
    end
    
    function [Dx Dy] = gradRound(X)        
        [M,N] = size(X);
        Dx = -diff(X,[],1);
        tmp = X(N,:)-X(1,:);
        Dx = [Dx;tmp];
        
        Dy = -diff(X,[],2);
        tmp = X(:,N)-X(:, 1);
        Dy = [Dy tmp];
    end
    