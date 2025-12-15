% The relationship between IED TW delays and XC delays，XC delays and XC strength,
% XC delays and Euclidean distance.
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
load(fullfile('step3_IEDtravelingWave', 'IEDtw_meth-None_TD-6d0.mat'));
% Load IED cross-correlation data with polynomial fitting
load(fullfile('step5_IEDxcorr_polyfit',  'IED6d0_xcorr100ms_Rmax_delay_allSubject.mat'));
% Process each subject individually
for nS = 1:length(subG)
        subID = subG{nS};
        % Load Euclidean distance matrix
        load(fullfile('step2_channelDistance',  [subID, '_EuclideabDistanceMatrix.mat']));
        %% Extract IED transfer matrix
        % Extract traveling wave delay matrix and transfer rate matrix
        IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
        IEDtwRate = IEDtw(nS).indirectTransferRateMatrix;
        % Extract cross-correlation delay and correlation coefficient matrices
        IEDxcD = Delay_allSubject{nS};
        IEDxcR = xcR_allSubject{nS};
        % Create a mask to identify valid data points and remove outliers
        % 1. Remove diagonal elements (self-comparisons)
        % 2. Exclude low transfer rates (<0.1)
        % 3. Exclude poor correlations (Fisher-z transformed threshold of 0.3)
        % 4. Exclude NaN/inf values in delays and correlations
        Mask = IEDtwRate;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask<0.1) = 0;
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(IEDxcR) | isinf(IEDxcR) | IEDxcR <= zThd) = 0;
        Mask(isnan(IEDxcD) | isinf(IEDxcD)) = 0;
        % Extract valid data points using the mask
        ED = DistanceWorld(Mask>0);
        IEDtwD = IEDtwD(Mask>0);
        IEDxcD = IEDxcD(Mask>0);
        IEDxcR = log(IEDxcR(Mask>0));
        %% Correlation analysis
        % Calculate Spearman correlation between traveling wave and cross-correlation delays
        [R_twxcD(nS), P_twxcD(nS)] = corr(IEDtwD, IEDxcD, "type", "Spearman");
        % Calculate Spearman correlation between cross-correlation delays and strength
        [R_xcDR(nS), P_xcDR(nS)] = corr(IEDxcD, IEDxcR, "type", "Spearman");
        % Analyze cross-correlation delays with polynomial fits (orders 1-6)
        for order = 1:6 
            % Fit polynomial model to remove distance effect from IED cross-correlation delays
            [polyModel, ~, Ftest] = f_polyfit(ED, IEDxcD, order);
            IEDxcDres = polyModel.residual;
            adjR2_xcDED(nS, order) = polyModel.adjR2;
            p_adjR2_xcDED(nS, order) = Ftest.p_value;
            % Remove distance effect from cross-correlation coefficients using linear fit
            [polyModel, ~, ~] = f_polyfit(ED, IEDxcR, order);
            IEDxcRres = polyModel.residual;
            % Calculate correlation between residuals of delays and coefficients
            [R_xcDR_rED(nS, order), P_xcDR_rED(nS, order)] = corr(IEDxcDres, IEDxcRres, "type", "Spearman");
            % Fit polynomial model to remove distance effect from IED traveling wave delays
            [polyModel, ~, ~] = f_polyfit(ED, IEDtwD, order);
            IEDtwDres = polyModel.residual;
            % Calculate correlation between residuals of traveling wave and cross-correlation delays
            [R_twxcD_rED(nS, order), P_twxcD_rED(nS, order)] = corr(IEDtwDres, IEDxcDres, "type", "Spearman");
        end
end
%% Save all analysis results
saveFolder = 'step5_IEDxcorr_polyfit';
saveName = 'IED6d0_XC100ms_polyfit_delays_Rmax_ED_es.mat';
save(fullfile(saveFolder, saveName), ...
    'R_twxcD', 'P_twxcD', 'R_twxcD_rED', 'P_twxcD_rED', ...
    'R_xcDR', 'P_xcDR', 'R_xcDR_rED', 'P_xcDR_rED', ...
    'adjR2_xcDED', 'p_adjR2_xcDED');