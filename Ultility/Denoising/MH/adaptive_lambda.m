function lambda = adaptive_lambda(cur_proj, A)
%
% A heuristic method of choosing an adaptive tikhonov regularization
% parameter.
%

std_Y = std(cur_proj);
std_A = std(A);
lam1 = std_Y^2/sum((std_A-std_Y).^2);
lam2 = lam1/std(std_A);
lambda = (lam1 + lam2)/2;