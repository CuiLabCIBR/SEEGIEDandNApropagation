% Perform bootstrapped Spearman correlation and polynomial regression 
% analyses to investigate relationships between IED traveling wave delays,
% non-IED cross-correlation delays, their correlation strengths, and 
% Euclidean distances, while controlling for distance effects through 
% parallel processing.
clc; clear; close all; addpath z_toolbox
load(fullfile('step6_NAxcorr_polyfit', ...
        'NAnon6d0IED_XC100ms_polyfit_delays_Rmas_ED_es.mat'), ...
        'EDas', 'IEDtwDas', 'NAxcDas', 'NAxcRas');
%% Parallel processing for bootstrap analysis (1000 iterations)
parfor nIf = 1:1000
        disp(nIf);
        % 1. Correlation between IED traveling wave delays and NA cross-correlation delays
        [R_twxcD_as(nIf), P_twxcD_as(nIf)] = iforest_spearman(IEDtwDas, NAxcDas);
        % 2. Same correlation after regressing out distance effect
        [R_twxcD_rED_as(nIf, :), P_twxcD_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(IEDtwDas, NAxcDas, EDas);
        % 3. Polynomial fit between cross-correlation delays and distance
        adjR2_xcDED_as(nIf, :) = iforest_polyfit(EDas, NAxcDas);
        % 4. Correlation between cross-correlation delays and strength
        [R_xcDR_as(nIf), P_xcDR_as(nIf)] = iforest_spearman(NAxcDas, NAxcRas);
        % 5. Same correlation after regressing out distance effect
        [R_xcDR_rED_as(nIf, :), P_xcDR_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(NAxcDas, NAxcRas, EDas);
end
%% Save all analysis results
saveFolder = 'step6_NAxcorr_polyfit';
saveName = 'NAnon6d0IED_XC100ms_polyfit_delays_Rmas_ED_as.mat';
save(fullfile(saveFolder, saveName), ...
    "R_twxcD_as", "P_twxcD_as", "R_twxcD_rED_as", "P_twxcD_rED_as", ...
    "adjR2_xcDED_as", ...
    "R_xcDR_as", "P_xcDR_as", "R_xcDR_rED_as", "P_xcDR_rED_as");