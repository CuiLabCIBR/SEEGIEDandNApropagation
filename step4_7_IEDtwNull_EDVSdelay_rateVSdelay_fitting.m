% Analyze null-model relationships among:
% 1) Euclidean distance (ED) and IED traveling-wave delay
% 2) IED traveling-wave transfer rate and delay
%
% For each null-model iteration, the script:
% - loads subject-level null-model traveling-wave matrices,
% - fits ED-delay relationships for each subject,
% - computes raw and ED-controlled correlations between rate and delay,
% - saves subject-level fitting and correlation results for all null iterations.

clc; clear; close all; 
addpath z_toolbox;
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

%% Load subject-level null-model data
% For each subject, load:
% 1) null-model traveling-wave delay/rate matrices
% 2) inter-channel Euclidean distance matrix
for nS = 1:length(subG)
    subID = subG{nS}; disp(subID);
    
    % Load subject-specific null-model traveling-wave results
    load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op3_IEDtwNull', [subID, '_IEDtw_mull_TD-6d0_IT-100ms.mat']));
    
    % Load Euclidean distance matrix for current subject
    load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
    distanceAll{nS} = DistanceWorld;

    % Organize null-model matrices by iteration
    for nNull = 1:1000
        delayAll{nNull}{nS} = IEDtwNull.indirectTransferDelayMatrix(:, :, nNull);
        rateAll{nNull}{nS} = IEDtwNull.indirectTransferRateMatrix(:, :, nNull);
    end
end

%% Fit null-model data for each permutation
parfor nNull = 1:1000
    disp(nNull);
    delayNull = delayAll{nNull};
    rateNull = rateAll{nNull};
    [ED_delay_expFit{nNull}, rate_delay_corr{nNull}] = subparfor(delayNull, rateNull, distanceAll);
end

%% Save results
saveFolder = 'step4_IEDtw_TD-6d0_IT-100ms_fitting';
mkdir(saveFolder);
saveName = 'op7_IEDtw_null_ED-vs-delay_rate-vs-delay_fitiing.mat';
save(fullfile(saveFolder, saveName),  'ED_delay_expFit', 'rate_delay_corr');

%% Local function for one null-model iteration
function [ED_delay_expFit, rate_delay_corr] = subparfor(delayNull, rateNull, distanceAll)
    for nS = 1:47
        IEDtwD = delayNull{nS};
        IEDtwR = rateNull{nS};
        DistanceWorld = distanceAll{nS};
        
        %% Build a mask to retain valid channel pairs
        % Exclude:
        % 1) diagonal elements
        % 2) invalid rate values (NaN)
        % 3) low transfer-rate pairs (< 0.1)
        Mask = IEDtwR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask<0.1) = 0;

        %% Extract valid observations
        ED = DistanceWorld(Mask>0);
        IEDtwD = IEDtwD(Mask>0);
        IEDtwR = IEDtwR(Mask>0);
        IEDtwR = log(IEDtwR);
        
        %% Correlation analysis between transfer rate and delay
        [R_rate_delay, P_rate_delay] = corr(IEDtwR, IEDtwD, "type", "Spearman");

        %% Fit delay as a function of Euclidean distance using an exponential model
        % Model: y = p1 + p2 * (1 - exp(-x / p3))
        modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x ./ p(3)));
        beta0 = [1; 1; 1];  alpha = 0.05;
        [~, res_ED_delay_exp, stats_ED_delay_exp, ~] = fit_lsqcurvefit(ED, IEDtwD, modelfun, beta0, alpha);

        %% Fit polynomial models and compute ED-controlled correlations
        R_rate_delay_rED_poly = [];
        P_rate_delay_rED_poly = [];
        R_rate_delay_rED_exp = [];
        P_rate_delay_rED_exp = [];

        for order = 1:10
            % Polynomial fit: ED -> delay
            [~, res_ED_delay_poly, ~, ~] = fit_polyfit(ED, IEDtwD, order, alpha);
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
        ED_delay_expFit{nS} = stats_ED_delay_exp;

        rate_delay_corr(nS).r = R_rate_delay;
        rate_delay_corr(nS).p = P_rate_delay;
        rate_delay_corr(nS).r_rED_exp = R_rate_delay_rED_exp;
        rate_delay_corr(nS).p_rED_exp = P_rate_delay_rED_exp;
        rate_delay_corr(nS).r_rED_poly = R_rate_delay_rED_poly;
        rate_delay_corr(nS).p_rED_poly = P_rate_delay_rED_poly;
    end
end