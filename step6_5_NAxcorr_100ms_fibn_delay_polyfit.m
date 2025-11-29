% Analyzes the relationship between fiber count and non-IED cross-correlation
% delays across multiple subjects using Spearman correlation, while 
% controlling for Euclidean distance effects through polynomial regression
% (orders 1-6) to remove distance confounding from both variables.
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
% Load fiber count data for all subjects
load(fullfile('step2_fiberTrack',  'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount', 'fiberLength');
% Load cross-correlation data
load(fullfile('step6_NAxcorr_polyfit',  'NAnon6d0IED_xcorr100ms_Rmax_delay_allSubject.mat'));
% Initialize arrays for aggregated data across all subjects
NAxcDas = []; % Aggregated cross-correlation delays
EDas = [];    % Aggregated Euclidean distances
fibNas = [];
% Process each subject individually
for nS = 1:length(subG)
        subID = subG{nS};
        % Load Euclidean distance matrix
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']));
        % Extract cross-correlation delay and correlation coefficient matrices
        NAxcD = Delay_allSubject{nS};
        NAxcR = xcR_allSubject{nS};
        fibN = fiberCount{nS};
        if isempty(fibN); continue; end
        Mask = NAxcR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(NAxcR) | isinf(NAxcR) | NAxcR <= zThd) = 0;
        Mask(isnan(NAxcD) | isinf(NAxcD)) = 0;
        Mask(fibN==0) = 0;
        % Extract valid data
        ED = DistanceWorld(Mask>0);
        NAxcD = NAxcD(Mask>0);
        fibN = fibN(Mask>0);
        fibN = log(fibN);
        % Aggregate data across subjects
        EDas = [EDas; ED];
        NAxcDas = [NAxcDas; NAxcD];
        fibNas = [fibNas; fibN];
        % Compute correlation between fiber count and cross-correlation delays
        [R_fibNxcD(nS), P_fibNxcD(nS)] = corr(fibN, NAxcD, "type", "Spearman");
        % Compute distance-residualized correlations
        for order = 1:6 
            % Remove distance effect from delays using polynomial fit
            [polyModel, ~, ~] = f_polyfit(ED, NAxcD, order);
            NAxcDres = polyModel.residual;
            % Remove distance effect from fiber counts using polynomial fit
            [polyModel, ~, ~] = f_polyfit(ED, fibN, order);
            fibNres = polyModel.residual;
            % Compute correlation between residuals
            [R_fibNxcD_rED(nS, order), P_fibNxcD_rED(nS, order)] = corr(fibNres, NAxcDres, "type", "Spearman");
        end
end
parfor nIf = 1:1000
        disp(nIf);
        % 1. Raw correlation between fiber count and IED cross-correlation
        [R_fibNxcD_as(nIf), P_fibNxcD_as(nIf)] = iforest_spearman(fibNas, NAxcDas);
        % 2. Distance-corrected correlation
        [R_fibNxcD_rED_as(nIf, :), P_fibNxcD_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(fibNas, NAxcDas, EDas);
end
% Save results
saveFolder = 'step6_NAxcorr_polyfit';
saveName = 'NAnon6d0IED_XC100ms_polyfit_delay_fibN.mat';
save(fullfile(saveFolder, saveName), ...
    'EDas', 'NAxcDas', 'fibNas', ...
    'R_fibNxcD', 'P_fibNxcD', 'R_fibNxcD_rED', 'P_fibNxcD_rED', ...
    'R_fibNxcD_as', 'P_fibNxcD_as', 'R_fibNxcD_rED_as', 'P_fibNxcD_rED_as');