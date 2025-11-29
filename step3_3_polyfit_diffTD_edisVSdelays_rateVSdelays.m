% Analyzes the relationship between IED propagation delays and rates while
% controlling for Euclidean distance effects, using Spearman correlation, 
% polynomial fitting, and bootstrap analysis for different IED detection thresholds.
clc;clear;close all; addpath z_toolbox;
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
% IED detection threshould
thd = 5:0.2:9;
thd_str = {'5d0', '5d2', '5d4', '5d6', '5d8', ...
        '6d0', '6d2', '6d4', '6d6', '6d8', ...
        '7d0', '7d2', '7d4', '7d6', '7d8', ...
        '8d0', '8d2', '8d4', '8d6', '8d8', '9d0'};
for nTD = 1:length(thd)
    % Load IED traveling wave data for current IED detection parameter
    load(fullfile('step3_IEDtravelingWave', ['IEDtw_meth-None_TD-', thd_str{nTD}, '.mat']));
    % Initialize aggregate variables for pooled subject data
    EDas = []; IEDtwDas = []; IEDtwRas = [];
    for nS = 1:length(subG)
        subID = subG{nS};
        % Load Euclidean distance matrix for subject
        load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
        % Extract IED traveling wave matrices
        IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
        IEDtwR = IEDtw(nS).indirectTransferRateMatrix;
        % Extract valid data points using mask
        Mask = IEDtwR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask<0.1) = 0;
        ED = DistanceWorld(Mask>0);
        IEDtwD = IEDtwD(Mask>0);
        IEDtwR = IEDtwR(Mask>0);
        IEDtwR = log(IEDtwR);
        %% Correlation Analysis
        % Raw correlation between log-rate and delay
        [R_twDR(nS), P_twDR(nS)] = corr(IEDtwR, IEDtwD, "type", "Spearman");
        for order = 1:6
            % Remove distance effect from log-rates using polynomial fit
            [polyModel, ~, ~] = f_polyfit(ED, IEDtwR, order);
            IEDtwR_rED = polyModel.residual;
            % Remove distance effect from delays using polynomial fit
            [polyModel, ~, ~] = f_polyfit(ED, IEDtwD, order);
            IEDtwD_rED = polyModel.residual;
            adjR2_twDED(nS, order) = polyModel.adjR2;
            % Correlation between residuals (rate vs delay)
            [R_twDR_rED(nS, order), P_twDR_rED(nS, order)] ...
                = corr(IEDtwR_rED, IEDtwD_rED, "type", "Spearman");
        end
        % Aggregate data across subjects
        EDas = [EDas; ED]; 
        IEDtwDas = [IEDtwDas; IEDtwD]; 
        IEDtwRas = [IEDtwRas; IEDtwR];
    end
    % Bootstrap analysis using parallel processing
    parfor nIf = 1:100
        disp([nTD, nIf]);
        % Bootstrap raw correlation between delays and rate
        [R_twDR_as(nIf), P_twDR_as(nIf)] = iforest_spearman(IEDtwDas, IEDtwRas);
        % Bootstrap residual correlations between delays and rate
        [R_twDR_rED_as(nIf, :), P_twDR_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(IEDtwDas, IEDtwRas, EDas);
        % Bootstrap polynomial fit quality between Euclidean distance and
        % delays
        adjR2_twDED_as(nIf, :) = iforest_polyfit(EDas, IEDtwDas);
    end
    % Save results for current IED detection paramater
    saveFolder = 'step3_IEDtw_polyfit_diffTD';
    mkdir(saveFolder);
    saveName = ['IEDtw_meth-none_TD-', thd_str{nTD}, '_polyfit.mat'];
    save(fullfile(saveFolder, saveName), ...
        'adjR2_twDED', ...
        'R_twDR', 'P_twDR', 'R_twDR_rED', 'P_twDR_rED', ...
        'adjR2_twDED_as', ...
        'R_twDR_as', 'P_twDR_as', 'R_twDR_rED_as', 'P_twDR_rED_as');
end