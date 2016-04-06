function [F, result] = KCS_GSR_V02(R, G, Y, opts, image1)
%       F = recSepTVNLM_Wx(R, G, Y, opts, image1)
%       function to recovery SepTV using Split Bregman method with seperable operator
%       under NLM constrainst
%       Input:
%           - R, G: sensing matrix size mxn, Phi*X = G'*X*R; or Phi =
%           Kornecker(G', R);
%           - Y: measurement matrix, mxm
%           - opts: parameter
%               + lambda: parameter for gradX and gradY
%               + mu: parameter for ||RFG - Y||_2
%               + nuy: parameter for V-->F
%               + gamma: Parameter for NLM
%               + nInner: maximum loop for Inner Loop, 15-25
%               + nOuter: number of maximum loop for Outer loop, 30-50
%               + tol: tolerance, 6e-7
%   Note:
%       - The total loop will be: nInner * nOuter, can be adjust between
%         these two parameters.
%       - NLM is performed in the end of each inner Loop, and was change to
%       boot up the result (0.05-.30 dB), by changing nInner, (When the
%       normal case does not gain), but it increase and then decrease-->
%       take care the "tol" value too.
% Implemented by Thuong Nguyen, 2013, April 03

%% test function
% clear all;
% load testLenna.mat;
%     nInner = 10;    nOuter = 60;       xSize = nInner;
%     tol = 5e-4;   gamma = 2;
%% Assign value;
% load('all_TV_lena_s0.1.mat');
nuy     = 0.05 ;
mu      = 1 ;
tol     = opts.par.tol;

% Reserve memory for the auxillary variables
Y0      = Y;
F       = opts.initial;

% Decompose the matrix R and G
[Ur, Lr] = eig(transpose(R)*R);         Ur = Ur';
[Ug, Lg] = eig(G*transpose(G));         Ug = Ug';

% The denomitor --> recheck
denoF = nuy + mu*diag(Lr,0)*diag(Lg,0)';

%% GSR based method
% load('all_TV_lena_s0.1.mat');
b = zeros(size(F));

iter = 0; 
err  = 1; 
display([ 'Init PSNR ' num2str(csnr(F, image1, 0, 0)) 'dB']);
while (iter < opts.par.IterNum && err > tol)
    % solving GSR sparse
    F_hat   = F;
    F_bar   = GSR_Solver_CS(F_hat - b, opts.par);
    
    % solving subproblem F    
    Fk      = Ur*(mu*(R)'*Y*(G)' + nuy*(F_bar + b))*Ug';    
    F       = (Ur)'*(Fk./denoF)*Ug;    
    Y       = Y+Y0-R*F*G;
    
    % update b
    b       = b - (F_hat - F_bar);
    
    iter        = iter + 1;
    err  = norm(F - F_hat,2)/norm(F_hat); 
    eachErr(iter) = err; 
    eachPSNR(iter) = csnr(F, image1, 0, 0); 
    eachSSIM(iter) = cal_ssim(F, image1, 0, 0); 
     
    % display info
    display(['iter out: ' num2str(iter) '/' num2str(opts.par.IterNum) ', PSNR:' ...
            num2str(csnr(F, image1, 0, 0)) 'dB, err:' num2str(err)]);
end

%% results
% result.eacherr = errorKPlus;
% if(opts.isShowPSNR == 1)
    result.eachPSNR = eachPSNR;
    result.eachErr = eachErr;
    %result.eachSSIM = eachSSIM;
% end



