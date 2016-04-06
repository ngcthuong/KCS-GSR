function [result] = AWTVNL1(Y, lambda, nuy, mu, nLoop, Wx, Wy, nl_spatial,par)
% [F, result] = ATV(Y, mu, nLoop, par)
%   Function for anisotropic SepTV recovery
%   - Input: 
%       + Y: Current measurement
%       + mu, lambda, nuy: are parameters of TV optimization problem
%       + nLoop: number of loop
%       + par: contain remainder parameters: R, dx, dy ....
%       + nl_spatial: the NLM filtered image in the spatial domain of the
%       previous iteration
%   - Output: 
%       + F: output image
%       + results: Contains other parameter: V, x, y, ... 


%% Parameters    
R       = par.R;        G       = par.G;      [M, N]    = size(R);      gamma   = par.gamma;
Ur      = par.Ur;       Lr      = par.Lr;     Ug        = par.Ug;       Lg      = par.Lg;
denoF   = par.denoF;    denoV   = par.denoV;

%% Initial value
if(par.init ==0)
    V   = par.V;        x       = par.x;        y       = par.y;
    bx  = par.bx;       by      = par.by;       W       = par.W;    
    F   = par.F;
else
    V   = zeros(N,N);   x       = V;            y       = V;
    bx  = V;            by      = V;            W       = V;
    F   = V;
end;
[dx dy] =gradCal3(F,3); dxt     = dx';          dyt     = dy';
[Ux,Lx] = eig(dxt*dx);  Ux      = Ux';
[Uy,Ly] = eig(dy*dyt);  Uy      = Uy';

%% Main loop
Yk      = mu*(R)'*Y*(G)';       
inner   = 0;  %F = zeros(N, N);

while(inner ~= nLoop)
    inner   = inner + 1;          F1 = F;
    
    % update F
    Fk      = Ur*(Yk + nuy*(V + W))*Ug';
    Fk      = Fk./denoF;
    F       = (Ur)'*Fk*Ug;
    
    % Update V
    
    if(inner == 1)
        Vk  = (Ux*(lambda*dxt*(x-bx) + lambda*(y - by)*dyt + nuy*F ...
                - nuy*W + gamma*nl_spatial )*transpose(Uy));
        Vk  = Vk./(denoV+gamma);
    else
        Vk  = (Ux*(lambda*dxt*(x-bx) + lambda*(y - by)*dyt + nuy*F ...
                    - nuy*W )*transpose(Uy));
        Vk  = Vk./(denoV);
    end    
    V       = (Ux)'*Vk*Uy;
    
    % update x and y
    grVx    = dx*(V) + bx;
    grVy    = (V)*dy + by;
    x       = max(abs(grVx) - Wx./lambda, 0).*sign(grVx);
    y       = max(abs(grVy) - Wy./lambda, 0).*sign(grVy);
    
    bx      = grVx-x;           by  = grVy-y;
    W       = V + W - F;
    
    % compare results;
    errK    = norm(F - F1,2)/norm(F,2);
    errAvg(inner) = errK;
%     figure(2);
%     imshow(F, []);
end

result.F = F;   result.V = V;   result.x = x;   result.y        = y;
result.bx=bx;   result.by=by;   result.W = W;   result.errAvg   = errAvg;

