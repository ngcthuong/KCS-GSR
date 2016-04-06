function [R, G] = KCS_SensingMtx(img_size, subrate, trial)
% function to generate sensing matrix
patch = ['SensingMatrixSepTV' num2str(img_size)];
if ~exist(patch, 'dir');
    display(['Generate new folder for sensing matrix']);
    mkdir(patch);
end;

projection_matrix_file 	= [patch '\SenMtr_R_trial' num2str(trial) '.mat'];
R                       = BCS_SPL_GenerateProjection_2(img_size, subrate, projection_matrix_file);

projection_matrix_file  = [patch '\SenMtr_G_trial' num2str(trial) '.mat'];
G                       = BCS_SPL_GenerateProjection_2(img_size, subrate, projection_matrix_file);
G                       = transpose(G);

