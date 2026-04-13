% Analyze the relationship between fiber count and IED traveling-wave delay
% for each subject.
%
% For each subject, the script:
% 1) extracts valid channel pairs,
% 2) computes the raw correlation between fiber count and delay,
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
% Load subject-level fiber-count matrices
load(fullfile('step2_fiberTrack', 'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
% Load subject-level IED traveling-wave results
load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op2_IEDtw_TD-6d0_IT-100ms.mat'));

for nS = 1:length(subG)
    subID = subG{nS}; disp(subID);
    % Load the inter-channel Euclidean distance matrix for the current subject
    load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
    %% Extract subject-specific fiber-count matrix
    fibN = fiberCount{nS};
    if isempty(fibN)
        continue;
    end
    
    %% Extract subject-specific traveling-wave matrices
    IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
    IEDtwR = IEDtw(nS).indirectTransferRateMatrix;

    %% Build a mask to retain valid channel pairs
    % Exclude:
    % 1) diagonal elements
    % 2) invalid rate values (NaN)
    % 3) low transfer-rate pairs (< 0.1)
    % 4) channel pairs with zero fiber count
    Mask = IEDtwR;
    Mask = Mask - diag(diag(Mask));
    Mask(isnan(Mask)) = 0;
    Mask(Mask<0.1) = 0;
    Mask(fibN==0) = 0;
    ED = DistanceWorld(Mask>0);
    IEDtwD = IEDtwD(Mask>0);
    fibN = fibN(Mask>0);
    fibN = log(fibN);

    %% Raw correlation between fiber count and traveling-wave delay
    [R_fibN_delay, P_fibN_delay] = corr(fibN, IEDtwD, "type", "Spearman");

   %% Fit delay as a function of Euclidean distance using an exponential model
    x = ED; y = IEDtwD; 
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x ./ p(3)));
    p1 = median(y, 'omitnan');
    p2 = (max(y, [], 'omitnan') - min(y, [], 'omitnan'))/2;
    p3 = max(std(x, 'omitnan'), 0.5);
    beta0 = [p1; p2; p3];
    alpha = 0.05;
    [~, res_ED_delay_exp, stats, ~] = fit_lsqcurvefit(x, y, modelfun, beta0, alpha);
    disp(['Adjusted R2 = ', num2str(stats.R2adj)]);
    
    %% Fit polynomial models and compute ED-controlled correlations
    R_fibN_delay_rED_poly = [];
    P_fibN_delay_rED_poly = [];
    R_fibN_delay_rED_exp = [];
    P_fibN_delay_rED_exp = [];
    for order = 1:10
        % Polynomial fit: ED -> delay
        [~,  res_ED_delay_poly, stats,  ~] = fit_polyfit(ED, IEDtwD, order, alpha);
        disp(['Adjusted R2 = ', num2str(stats.R2adj)]);
        % Polynomial fit: ED -> fiber count
        [~, res_ED_fibN_poly, ~, ~] = fit_polyfit(ED, fibN, order, alpha);
        % Correlation between polynomial residuals of fiber count and delay
        [R_fibN_delay_rED_poly(order),  P_fibN_delay_rED_poly(order)] ...
            = corr(res_ED_fibN_poly, res_ED_delay_poly,  "type", "Spearman");
        % Correlation between polynomial residuals of fiber count and
        % exponential residuals of delay
        [R_fibN_delay_rED_exp(order), P_fibN_delay_rED_exp(order)] ...
            = corr(res_ED_fibN_poly, res_ED_delay_exp,  "type", "Spearman");
    end

    %% Save subject-level results
    fibN_delay_corr(nS).delay = IEDtwD;
    fibN_delay_corr(nS).fibN = fibN;
    fibN_delay_corr(nS).ED = ED;
    fibN_delay_corr(nS).r = R_fibN_delay;
    fibN_delay_corr(nS).p = P_fibN_delay;
    fibN_delay_corr(nS).r_rED_exp = R_fibN_delay_rED_exp;
    fibN_delay_corr(nS).p_rED_exp = P_fibN_delay_rED_exp;
    fibN_delay_corr(nS).r_rED_poly = R_fibN_delay_rED_poly;
    fibN_delay_corr(nS).p_rED_poly = P_fibN_delay_rED_poly;
end

%% Save all results
saveFolder = 'step4_IEDtw_TD-6d0_IT-100ms_fitting';
mkdir(saveFolder);
saveName = 'op6_IEDtw_fibN-VS-delay_fitting.mat';
save(fullfile(saveFolder, saveName), 'fibN_delay_corr');