% Analyze the relationships among:
% 1) Euclidean distance (ED) and IED traveling-wave delay
% 2) IED traveling-wave transfer rate and delay
%
% For each subject, the script:
% - extracts valid channel pairs,
% - fits ED-delay relationships using exponential and polynomial models,
% - computes the raw correlation between transfer rate and delay,
% - computes ED-controlled correlations using model residuals,
% - saves subject-level fitting and correlation results.
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
% Load subject-level IED traveling-wave results
load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting',  'op2_IEDtw_TD-6d0_IT-100ms.mat'));
for nS = 1:length(subG)
        subID = subG{nS}; disp(subID);
        % Load inter-channel Euclidean distance matrix for the current subject
        load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
        
        %% Extract subject-specific traveling-wave matrices
        % IEDtwD: traveling-wave delay matrix
        % IEDtwR: traveling-wave transfer-rate matrix
        IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
        IEDtwR = IEDtw(nS).indirectTransferRateMatrix;

        %% Build a mask to retain valid channel pairs
        % Exclude:
        % 1) diagonal elements (self-self pairs)
        % 2) invalid rate values (NaN)
        % 3) low transfer-rate pairs (< 0.1)
        Mask = IEDtwR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask<0.1) = 0;
        ED = DistanceWorld(Mask>0);

        %% Extract valid observations
        IEDtwD = IEDtwD(Mask>0);
        IEDtwR = IEDtwR(Mask>0);
        IEDtwR = log(IEDtwR);

        %% Correlation analysis between transfer rate and delay
        [R_rate_delay, P_rate_delay] = corr(IEDtwR, IEDtwD, "type", "Spearman");

        %% Fit delay as a function of Euclidean distance using an exponential model
        % Model: y = p1 + p2 * (1 - exp(-x / p3))
        x = ED; y = IEDtwD; 
        modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x ./ p(3)));
        p1 = median(y, 'omitnan');
        p2 = (max(y, [], 'omitnan') - min(y, [], 'omitnan'))/2;
        p3 = max(std(x, 'omitnan'), 0.5);
        beta0 = [p1; p2; p3];
        alpha = 0.05;
        [beta_ED_delay_exp, res_ED_delay_exp, stats_ED_delay_exp, curve_ED_delay_exp] ...
            = fit_lsqcurvefit(x, y, modelfun, beta0, alpha);
        disp(['Adjusted R2 = ', num2str(stats_ED_delay_exp.R2adj)]);
        %% Fit polynomial models and compute ED-controlled correlations
        beta_ED_delay_poly = {};
        stats_ED_delay_poly = {};
        curve_ED_delay_poly = {};

        R_rate_delay_rED_poly = [];
        P_rate_delay_rED_poly = [];
        R_rate_delay_rED_exp = [];
        P_rate_delay_rED_exp = [];

        for order = 1:10
            % Polynomial fit: ED -> delay
            [beta_ED_delay_poly{order}, res_ED_delay_poly, ...
                stats_ED_delay_poly{order}, curve_ED_delay_poly{order}] ...
                = fit_polyfit(ED, IEDtwD, order, alpha);
            disp(['Adjusted R2 = ', num2str(stats_ED_delay_poly{order}.R2adj)]);
            
            % Polynomial fit: ED -> transfer rate
            [~, res_ED_rate_poly, ~, ~] =  fit_polyfit(ED, IEDtwR, order, alpha);
            
            % Correlation between polynomial residuals of rate and delay
            [R_rate_delay_rED_poly(order), P_rate_delay_rED_poly(order)] = ...
                corr(res_ED_rate_poly, res_ED_delay_poly,  "type", "Spearman");
            
            % Correlation between polynomial residuals of rate and
            % exponential residuals of delay
            [R_rate_delay_rED_exp(order), P_rate_delay_rED_exp(order)] = ...
                corr(res_ED_rate_poly, res_ED_delay_exp,  "type", "Spearman");
        end
        
        %% Save subject-level outputs
        ED_delay_fit(nS).ED = ED;
        ED_delay_fit(nS).delay = IEDtwD;
        ED_delay_fit(nS).beta_exp = beta_ED_delay_exp;
        ED_delay_fit(nS).stats_exp = stats_ED_delay_exp;
        ED_delay_fit(nS).curve_exp = curve_ED_delay_exp;
        ED_delay_fit(nS).beta_poly = beta_ED_delay_poly;
        ED_delay_fit(nS).stats_poly = stats_ED_delay_poly;
        ED_delay_fit(nS).curve_poly = curve_ED_delay_poly;

        rate_delay_corr(nS).delay = IEDtwD;
        rate_delay_corr(nS).rate = IEDtwR;
        rate_delay_corr(nS).ED = ED;
        rate_delay_corr(nS).r = R_rate_delay;
        rate_delay_corr(nS).p = P_rate_delay;
        rate_delay_corr(nS).r_rED_exp = R_rate_delay_rED_exp;
        rate_delay_corr(nS).p_rED_exp = P_rate_delay_rED_exp;
        rate_delay_corr(nS).r_rED_poly = R_rate_delay_rED_poly;
        rate_delay_corr(nS).p_rED_poly = P_rate_delay_rED_poly;
end
%% Save all results
saveFolder = 'step4_IEDtw_TD-6d0_IT-100ms_fitting';
mkdir(saveFolder);
saveName = 'op5_IEDtw_ED-vs-delay_rate-vs-delay_fitting.mat';
save(fullfile(saveFolder, saveName), "ED_delay_fit", "rate_delay_corr");