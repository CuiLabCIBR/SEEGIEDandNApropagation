% Analyze the relationships between IED traveling-wave delays, non-IED
% cross-correlation delays, cross-correlation strengths, and Euclidean
% distance for each subject using Spearman correlation and polynomial
% regression (orders 1-10) to control for distance effects.
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
%% Load subject-level traveling-wave results
 load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op2_IEDtw_TD-6d0_IT-100ms.mat'));
%% Load subject-level non-IED cross-correlation results
load(fullfile('step6_NAxcorr_TD-6d0_IT-100ms_fitting', 'op1_NAxcorr_delay_FC.mat'));
%% Process each subject individually
for nS = 1:length(subG)
        subID = subG{nS};
        % Load Euclidean distance matrix
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']));
        %% Extract subject-specific matrices
        % Traveling-wave delay and transfer-rate matrices
        twD = IEDtw(nS).indirectTransferDelayMatrix;
        twRate = IEDtw(nS).indirectTransferRateMatrix;
        % Non-IED cross-correlation delay and strength matrices
        xcD = NAxcDelay_as{nS};
        xcFC = NAxcFC_as{nS};

        %% Build a mask to retain valid channel pairs
        % Exclude:
        % 1) diagonal elements
        % 2) low traveling-wave transfer rate (< 0.1)
        % 3) weak XC strength values (threshold corresponding to r = 0.3 in Fisher-z space)
        % 4) invalid values in XC delay / XC strength / TW delay
        Mask = twRate;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask<0.1) = 0;
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(xcFC) | isinf(xcFC) | xcFC <= zThd) = 0;
        Mask(isnan(xcD) | isinf(xcD)) = 0;
        %% Extract valid observations
        ED = DistanceWorld(Mask>0);
        twD = twD(Mask>0);
        xcD = xcD(Mask>0);
        xcFC = log(xcFC(Mask>0));
        
        %% Fit XC delay as a function of Euclidean distance
        % y = p1 + p2 * (1 - exp(-x / p3))
        modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x ./ p(3)));
        beta0 = [1; 1; 1]; 
        alpha = 0.05;
        [beta_ED_xcD_exp, res_ED_xcD_exp,  ...
            stats_ED_xcD_exp, curve_ED_xcD_exp] ...
            = fit_lsqcurvefit(ED, xcD, modelfun, beta0, alpha);
        
        % Polynomial fits (orders 1 to 10)
        beta_ED_xcD_poly = {}; 
        res_ED_xcD_poly = [];
        stats_ED_xcD_poly = {}; 
        curve_ED_xcD_poly = {};
        for order = 1:10
            [beta_ED_xcD_poly{order}, res_ED_xcD_poly(:, order), ...
                stats_ED_xcD_poly{order}, curve_ED_xcD_poly{order}] ...
                = fit_polyfit(ED, xcD, order, alpha);
        end

        % Store fitting results
        ED_xcD_fit(nS).ED = ED;
        ED_xcD_fit(nS).xcD = xcD;
        ED_xcD_fit(nS).beta_exp = beta_ED_xcD_exp;
        ED_xcD_fit(nS).stats_exp = stats_ED_xcD_exp;
        ED_xcD_fit(nS).curve_exp = curve_ED_xcD_exp;
        ED_xcD_fit(nS).beta_poly = beta_ED_xcD_poly;
        ED_xcD_fit(nS).stats_poly = stats_ED_xcD_poly;
        ED_xcD_fit(nS).curve_poly = curve_ED_xcD_poly;

        %% Correlation between XC delay and XC strength
        [R_xcFC_xcD, P_xcFC_xcD] = corr(xcD, xcFC, "type", "Spearman");

        R_xcFC_xcD_rED_poly = []; 
        P_xcFC_xcD_rED_poly = [];
        R_xcFC_xcD_rED_exp = []; 
        P_xcFC_xcD_rED_exp = [];
        for order = 1:10
            [~, res_ED_xcFC_poly, ~, ~] =  fit_polyfit(ED, xcFC, order, alpha);
            % Correlation between polynomial residuals of XC strength and XC delay
            [R_xcFC_xcD_rED_poly(order), ...
                P_xcFC_xcD_rED_poly(order)] = ...
                corr(res_ED_xcFC_poly, res_ED_xcD_poly(:, order),  "type", "Spearman");
            % Correlation between polynomial residuals of XC strength and
            % exponential residuals of XC delay
            [R_xcFC_xcD_rED_exp(order), ...
                P_xcFC_xcD_rED_exp(order)] = ...
                corr(res_ED_xcFC_poly, res_ED_xcD_exp,  "type", "Spearman");
        end

        % Store correaltion results
        xcFC_xcD_corr(nS).xcD = xcD;
        xcFC_xcD_corr(nS).xcFC = xcFC;
        xcFC_xcD_corr(nS).r = R_xcFC_xcD;
        xcFC_xcD_corr(nS).p = P_xcFC_xcD;
        xcFC_xcD_corr(nS).r_rED_exp = R_xcFC_xcD_rED_exp;
        xcFC_xcD_corr(nS).p_rED_exp = P_xcFC_xcD_rED_exp;
        xcFC_xcD_corr(nS).r_rED_poly = R_xcFC_xcD_rED_poly;
        xcFC_xcD_corr(nS).p_rED_poly = P_xcFC_xcD_rED_poly;

        %% Correlation between XC delay and TW delay
        % Direct Spearman correlation
        [R_xcD_twD, P_xcD_twD] = corr(xcD, twD, "type", "Spearman");
        
        % Remove ED-related effects using exponential fitting
        [~, res_ED_twD_exp, ~, ~] ...
            = fit_lsqcurvefit(ED, twD, modelfun, beta0, alpha);

        [R_xcD_twD_rED_exp, P_xcD_twD_rED_exp] ...
            = corr(res_ED_xcD_exp, res_ED_twD_exp, "type", "Spearman");
        
        % Remove ED-related effects using polynomial fitting
        R_xcD_twD_rED_poly = []; 
        P_xcD_twD_rED_poly = [];
        for order =1:10
            [~, res_ED_twD_poly, ~, ~] = fit_polyfit(ED, twD, order, alpha);
            
            [R_xcD_twD_rED_poly(order), P_xcD_twD_rED_poly(order)] = ...
                corr(res_ED_xcD_poly(:, order), res_ED_twD_poly,  "type", "Spearman");
        end

        % Store correaltion results
        xcD_twD_corr(nS).xcD = xcD;
        xcD_twD_corr(nS).twD = twD;
        xcD_twD_corr(nS).r = R_xcD_twD;
        xcD_twD_corr(nS).p = P_xcD_twD;
        xcD_twD_corr(nS).r_rED_exp = R_xcD_twD_rED_exp;
        xcD_twD_corr(nS).p_rED_exp = P_xcD_twD_rED_exp;
        xcD_twD_corr(nS).r_rED_poly = R_xcD_twD_rED_poly;
        xcD_twD_corr(nS).p_rED_poly = P_xcD_twD_rED_poly;
end
%% Save all analysis results
saveFolder = 'step6_NAxcorr_TD-6d0_IT-100ms_fitting';
saveName = 'op2_NAxcorr_ED-vs-xcD_xcFC-vs-xcD.mat';
save(fullfile(saveFolder, saveName), 'ED_xcD_fit', 'xcFC_xcD_corr', 'xcD_twD_corr');