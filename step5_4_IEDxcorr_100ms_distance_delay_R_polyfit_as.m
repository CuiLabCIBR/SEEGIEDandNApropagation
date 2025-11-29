% Perform bootstrapped Spearman correlation and polynomial regression 
% analyses to assess relationships between IED traveling wave delays, 
% cross-correlation delays, their strengths, and distance effects. 
clc;clear;close all; addpath z_toolbox
load(fullfile('step5_IEDxcorr_polyfit', ...
                    'IED6d0_XC100ms_polyfit_delays_Rmax_ED_es.mat'), ...
                    'EDas', 'IEDtwDas', 'IEDxcDas', 'IEDxcRas');
%% Parallel processing for bootstrap analysis (1000 iterations)
parfor nIf = 1:1000
        disp(nIf);
        % 1. Correlation between IED traveling wave delays and cross-correlation delays
        [R_twxcD_as(nIf), P_twxcD_as(nIf)] = iforest_spearman(IEDtwDas, IEDxcDas);
        % 2. Same correlation after regressing out distance effect
        [R_twxcD_rED_as(nIf, :), P_twxcD_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(IEDtwDas, IEDxcDas, EDas);
        % 3. Polynomial fit between cross-correlation delays and distance
        adjR2_xcDED_as(nIf, :) = iforest_polyfit(EDas, IEDxcDas);
        % 4. Correlation between cross-correlation delays and strength
        [R_xcDR_as(nIf), P_xcDR_as(nIf)] = iforest_spearman(IEDxcRas, IEDxcDas);
        % 5. Same correlation after regressing out distance effect
        [R_xcDR_rED_as(nIf, :), P_xcDR_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(IEDxcDas, IEDxcRas, EDas);
end
%% Save all analysis results
saveFolder = 'step5_IEDxcorr_polyfit';
saveName = 'IED6d0_XC100ms_polyfit_delay_Rmax_ED_as.mat';
save(fullfile(saveFolder, saveName), ...
    'R_twxcD_as', 'P_twxcD_as', 'R_twxcD_rED_as', 'P_twxcD_rED_as', ...
    'adjR2_xcDED_as', ...
    'R_xcDR_as', 'P_xcDR_as', 'R_xcDR_rED_as', 'P_xcDR_rED_as');