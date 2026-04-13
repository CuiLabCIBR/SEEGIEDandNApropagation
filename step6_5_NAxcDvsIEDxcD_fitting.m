% Analyze the relationship between non-IED cross-correlation delay (NAxcD)
% and IED cross-correlation delay (IEDxcD) for each subject.
%
% For each subject, the script:
% 1) extracts valid channel pairs shared by NAxcorr and IEDxcorr,
% 2) computes the raw Spearman correlation between NAxcD and IEDxcD,
% 3) controls for Euclidean distance (ED) using exponential and polynomial residuals,
% 4) saves subject-level results for downstream analysis.
clc; clear; close all; addpath z_toolbox;
subG = {'sub-0001', 'sub-0002', 'sub-0003', 'sub-0004', 'sub-0005', ...
            'sub-0007', 'sub-0008', 'sub-0009', 'sub-0010', 'sub-0011', ...
            'sub-0012', 'sub-0013', 'sub-0014', 'sub-0015', 'sub-0016', ...
            'sub-0017', 'sub-0019', 'sub-0020', 'sub-0021', 'sub-0022', ...
            'sub-0023', 'sub-0028', 'sub-0029', 'sub-0031', 'sub-0033', ...
            'sub-0034', 'sub-0035', 'sub-0036', 'sub-0039', 'sub-0046', ...
            'sub-0047', 'sub-0050', 'sub-0052', 'sub-0054', 'sub-0062', ...
            'sub-0063', 'sub-0066', 'sub-0067', 'sub-0074', 'sub-0077', ...
            'sub-0087', 'sub-0089', 'sub-0092', 'sub-0093', 'sub-0094', ...
            'sub-0113', 'sub-0115'};
%% Load subject-level cross-correlation results
load(fullfile('step5_IEDxcorr_TD-6d0_IT-100ms_fitting',  'op1_IEDxcorr_delay_FC_allSubject.mat'));
%% Load subject-level non-IED cross-correlation results
load(fullfile('step6_NAxcorr_TD-6d0_IT-100ms_fitting', 'op1_NAxcorr_delay_FC.mat'));
%% Process each subject individually
for nS = 1:length(subG)
    subID = subG{nS};
    % Load Euclidean distance matrix for the current subject
    load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']));
    %% Extract subject-specific matrices
    IEDxcD = IEDxcDelay_as{nS};
    IEDxcFC = IEDxcFC_as{nS};
    NAxcD = NAxcDelay_as{nS};
    NAxcFC = NAxcFC_as{nS};

    %% Build a mask to retain valid channel pairs
    % Exclude:
    % 1) diagonal elements
    % 2) weak IED or non-IED XC strength values
    % 3) invalid delay values
    Mask = IEDxcFC;
    Mask = Mask - diag(diag(Mask));
    zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
    Mask(isnan(NAxcFC) | isinf(NAxcFC) | NAxcFC <= zThd) = 0;
    Mask(isnan(NAxcD) | isinf(NAxcD)) = 0;
    Mask(isnan(IEDxcFC) | isinf(IEDxcFC) | IEDxcFC <= zThd) = 0;
    Mask(isnan(IEDxcD) | isinf(IEDxcD)) = 0;
    %% Extract valid observations
    ED = DistanceWorld(Mask>0);
    IEDxcD = IEDxcD(Mask>0);
    NAxcD = NAxcD(Mask>0);
    
    %% Raw correlation between NAxcD and IEDxcD
    [R_NAxcD_IEDxcD, P_NAxcD_IEDxcD] = corr(NAxcD, IEDxcD, 'Type', 'Spearman');
    %% Fit XC delay as a function of Euclidean distance
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x ./ p(3)));
    beta0 = [1; 1; 1]; 
    alpha = 0.05;
    [~, res_ED_IEDxcD_exp, ~, ~] = fit_lsqcurvefit(ED, IEDxcD, modelfun, beta0, alpha);
    [~, res_ED_NAxcD_exp, ~, ~] = fit_lsqcurvefit(ED, NAxcD, modelfun, beta0, alpha);
    [R_NAxcD_IEDxcD_rED_exp, P_NAxcD_IEDxcD_rED_exp] = ...
            corr(res_ED_NAxcD_exp, res_ED_IEDxcD_exp,  "type", "Spearman");
    
    %% Polynomial residual correlations after removing ED effects
    R_xcFC_xcD_rED_poly = []; 
    P_xcFC_xcD_rED_poly = [];
    for order = 1:10
        [~, res_ED_NAxcD_poly, ~, ~] =  fit_polyfit(ED, NAxcD, order, alpha);
        [~, res_ED_IEDxcD_poly, ~, ~] =  fit_polyfit(ED, IEDxcD, order, alpha);
        % Correlation between polynomial residuals of XC strength and XC delay
        [R_NAxcD_IEDxcD_rED_poly(order), P_NAxcD_IEDxcD_rED_poly(order)] = ...
            corr(res_ED_NAxcD_poly, res_ED_IEDxcD_poly,  "type", "Spearman");
    end

    %% Save subject-level outputs
    NAxcD_IEDxcD_corr(nS).NAxcD = NAxcD;
    NAxcD_IEDxcD_corr(nS).IEDxcD = IEDxcD;
    NAxcD_IEDxcD_corr(nS).r = R_NAxcD_IEDxcD;
    NAxcD_IEDxcD_corr(nS).p = P_NAxcD_IEDxcD;
    NAxcD_IEDxcD_corr(nS).r_rED_exp = R_NAxcD_IEDxcD_rED_exp;
    NAxcD_IEDxcD_corr(nS).p_rED_exp = P_NAxcD_IEDxcD_rED_exp;
    NAxcD_IEDxcD_corr(nS).r_rED_poly = R_NAxcD_IEDxcD_rED_poly;
    NAxcD_IEDxcD_corr(nS).p_rED_poly = P_NAxcD_IEDxcD_rED_poly;
end
%% Save all analysis results
saveFolder = 'step6_NAxcorr_TD-6d0_IT-100ms_fitting';
saveName = 'op4_NAxcorr_NAxcD-vs-IEDxcD.mat';
save(fullfile(saveFolder, saveName), 'NAxcD_IEDxcD_corr');