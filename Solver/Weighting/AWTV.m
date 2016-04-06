function [result] = AWTV(Y, lambda, nuy, mu, nLoop, Wx, Wy, par)
% [F, result] = ATV(Y, mu, nLoop, par)
%   Function for anisotropic SepTV recovery
%   - Input: 
%       + Y: Current measurement
%       + mu, lambda, nuy: are parameters of TV optimization problem
%       + nLoop: number of loop
%       + par: contain remainder parameters: R, dx, dy ....
%   - Output: 
%       + F: output image
%       + results: Contains other parameter: V, x, y, ... 


%% Parameters    
R = par.R;      G = par.G;      N = max(size(R));
dx = par.dx;    dy = par.dy;    dxt = dx';      dyt = dy';
Ur = par.Ur;    Ug = par.Ug;    
% Ux = par.Ux;    Uy = par.Uy;    
denoF = par.denoF;
denoV = par.denoV;
%% Initial value
if(par.init ==0)
    V = par.V;    x = par.x;    y = par.y;
    bx= par.bx;   by = par.by;   W = par.W;    
    F = par.F;
else
    V = zeros(N,N);  x= V;       y = V;
    bx = V;         by = V;     W = V;
    F = V;
end;
% [dx dy] = gradCal3(F, 3);
%  dxt = dx';      dyt = dy';
%  [Ur, Lr] = eig(transpose(R)*R);         Ur = Ur';
%  [Ug, Lg] = eig(G*transpose(G));         Ug = Ug';
 [Ux, Lx] = eig(dxt*dx);                 Ux = Ux';
 [Uy, Ly] = eig(dy*dyt);                 Uy = Uy';
%  denoF = nuy + mu*diag(Lr,0)*diag(Lg,0)';
%  denoV = nuy + lambda* (diag(Lx,0)*ones(1,N) + (diag(Ly,0)*ones(1,N))' );

%% Main loop
Yk = mu*(R)'*Y*(G)';
inner = 0;  %F = zeros(N, N);

while(inner ~= nLoop)
    inner = inner + 1;          %F1 = F;
    
    % update F
    Fk = Ur*(Yk + nuy*(V + W))*Ug';
    Fk = Fk./denoF;
    F = (Ur)'*Fk*Ug;
    
    % Update V
    Vk = (Ux*(lambda* (dxt*(x-bx) + (y - by)*dyt) + nuy*(F - W))*(Uy)');
    Vk = Vk./denoV;
    V = (Ux)'*Vk*Uy;
    
    % update x and y
    grVx = dx*(V) + bx;
    grVy = (V)*dy + by;
    x = max(abs(grVx) - Wx./lambda, 0).*sign(grVx);
    y = max(abs(grVy) - Wy./lambda, 0).*sign(grVy);
    
    bx = grVx-x;
    by = grVy-y;
    W = V + W - F;
    
    % compare results;
%     errorKPlus = norm(F - F1,2)/norm(F,2);
%     errAvg(inner) = errorKPlus;
%     figure(2);
%     imshow(F, []);
end
%     Y = Y+Y0-R*F*G;
result.F = F;   result.V = V;
result.x = x;   result.y = y;
result.bx=bx;   result.by=by;   result.W = W;
% result.errAvg = errAvg;

