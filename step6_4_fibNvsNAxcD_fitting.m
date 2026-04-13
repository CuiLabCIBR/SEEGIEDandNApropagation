% Analyze the relationship between fiber count and non-IED cross-correlation
% delay for each subject, while controlling for the effect of
% inter-channel Euclidean distance (ED) using residuals from polynomial
% or exponential regression models.
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
%% Load subject-level non-IED cross-correlation results
load(fullfile('step6_NAxcorr_TD-6d0_IT-100ms_fitting', 'op1_NAxcorr_delay_FC.mat'));
%% Load subject-level structural connectivity results (fiber count)
load(fullfile('step2_fiberTrack',  'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
%% Process each subject individually
fibNas = []; xcDas = []; xcDres_as = []; fibNres_as = [];
for nS = 1:length(subG)
        subID = subG{nS};
        % Load Euclidean distance matrix
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']));
        
        %% Extract subject-specific matrices
        xcD = NAxcDelay_as{nS};
        xcR = NAxcFC_as{nS};
        fibN = fiberCount{nS};
        if isempty(fibN)
            continue;
        end
        
        %% Build a mask to retain valid channel pairs
        % Exclude:
        % 1) diagonal elements
        % 2) weak XC strength values
        % 3) invalid XC delay values
        % 4) channel pairs with zero fiber count
        Mask = xcR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(xcR) | isinf(xcR) | xcR <= zThd) = 0;
        Mask(isnan(xcD) | isinf(xcD)) = 0;
        Mask(fibN==0) = 0;
        %% Extract valid observations
        ED = DistanceWorld(Mask>0);
        xcD = xcD(Mask>0);
        fibN = fibN(Mask>0);
        fibN = log(fibN);

        %% Raw correlation between fiber count and XC delay
        [R_fibN_xcD, P_fibN_xcD] = corr(fibN, xcD, "type", "Spearman");
        %% Fit XC delay as a function of ED using an exponential model
        modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x ./ p(3)));
        beta0 = [1; 1; 1]; 
        alpha = 0.05;
        [~, res_ED_xcD_exp, ~, ~] = fit_lsqcurvefit(ED, xcD, modelfun, beta0, alpha);
    
        %% Fit polynomial models to remove ED-related effects
        R_fibN_xcD_rED_poly = []; P_fibN_xcD_rED_poly = [];
        R_fibN_xcD_rED_exp = []; P_fibN_xcD_rED_exp = [];
        res_ED_fibN_poly = [];
        for order = 1:10
            % Polynomial fit: ED -> XC delay
            [~,  res_ED_xcD_poly, ~] = fit_polyfit(ED, xcD, order, alpha);
            % Polynomial fit: ED -> fiber count
            [~, res_ED_fibN_poly(:, order), ~, ~] = fit_polyfit(ED, fibN, order, alpha);
            % Correlation between polynomial residuals
            [R_fibN_xcD_rED_poly(order),  P_fibN_xcD_rED_poly(order)] ...
                = corr(res_ED_fibN_poly(:, order), res_ED_xcD_poly,  "type", "Spearman");
            % Correlation between fiber-count polynomial residuals and
            % XC-delay exponential residuals
            [R_fibN_xcD_rED_exp(order), P_fibN_xcD_rED_exp(order)] ...
                = corr(res_ED_fibN_poly(:, order), res_ED_xcD_exp,  "type", "Spearman");
        end

        %% Save outputs
        fibN_xcD_corr(nS).xcD = xcD;
        fibN_xcD_corr(nS).fibN = fibN;
        fibN_xcD_corr(nS).r = R_fibN_xcD;
        fibN_xcD_corr(nS).p = P_fibN_xcD;
        fibN_xcD_corr(nS).r_rED_exp = R_fibN_xcD_rED_exp;
        fibN_xcD_corr(nS).p_rED_exp = P_fibN_xcD_rED_exp;
        fibN_xcD_corr(nS).r_rED_poly = R_fibN_xcD_rED_poly;
        fibN_xcD_corr(nS).p_rED_poly = P_fibN_xcD_rED_poly;
        %%
        fibNas = [fibNas; fibN];  xcDas = [xcDas; xcD];
        xcDres_as = [xcDres_as; res_ED_xcD_exp];
        fibNres_as = [fibNres_as; res_ED_fibN_poly(:, 6)];
end
% Save results
saveFolder = 'step6_NAxcorr_TD-6d0_IT-100ms_fitting';
saveName = 'op3_NAxcorr_fibN-vs-xcD_fitting.mat';
save(fullfile(saveFolder, saveName), 'fibN_xcD_corr', ...
    'fibNas', 'xcDas', 'xcDres_as', 'fibNres_as');