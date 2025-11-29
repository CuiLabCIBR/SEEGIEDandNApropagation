% Analyzes the relationship between white matter fiber count and IED
% propagation delays while controlling for spatial distance effects 
% through polynomial residualization and bootstrap validation.
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
% Load fiber matrix data containing fiber count information
load(fullfile('step2_fiberTrack', 'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
TD = 5:0.2:9;
TDstr = {'5d0', '5d2', '5d4', '5d6', '5d8', ...
                '6d0', '6d2', '6d4', '6d6', '6d8', ...
                '7d0', '7d2', '7d4', '7d6', '7d8', ...
                '8d0', '8d2', '8d4', '8d6', '8d8', '9d0'};
for nTD = 1:length(TD)
    % Load IED traveling wave data for current IED detection parameter
    load(fullfile('step8_IEDtravelingWave_corrSelect', ['IEDtw_meth-corrSelect_TD-', TDstr{nTD}, '.mat']));
    EDas = []; IEDtwDas = []; fibNas = [];
    for nS = 1:length(subG)
        subID = subG{nS};
        load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
        fibN = fiberCount{nS};
        if isempty(fibN)
            continue;
        end
        IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
        IEDtwR = IEDtw(nS).indirectTransferRateMatrix;
        % Extract relevant values using mask
        Mask = IEDtwR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(fibN==0) = 0;
        ED = DistanceWorld(Mask>0);
        IEDtwD = IEDtwD(Mask>0);
        fibN = fibN(Mask>0);
        fibN = log(fibN);
        %% Correlation Analysis: fiber count vs. IED transfer delay
        [R_twDfibN(nS), P_twDfibN(nS)] = corr(fibN, IEDtwD, "type", "Spearman");
        for order = 1:6
            % Remove Euclidean distance effect from fiber count using polynomial fit
            [polyModel, ~, ~] = f_polyfit(ED, fibN, order);
            fibNrED = polyModel.residual;
            % Remove Euclidean distance effect from IED delay using polynomial fit
            [polyModel, ~, ~] = f_polyfit(ED, IEDtwD, order);
            IEDtwDrED = polyModel.residual;
            % Calculate correlation between residuals
            [R_twDfibN_rED(nS, order), P_twDfibN_rED(nS, order)] = corr(fibNrED, IEDtwDrED, "type", "Spearman");
        end
        % Aggregate data across subjects for bootstrap analysis
        EDas = [EDas; ED]; 
        IEDtwDas = [IEDtwDas; IEDtwD];
        fibNas = [fibNas; fibN];
    end
    %% Bootstrap analysis using parallel processing
    parfor nIf = 1:100
        disp([nTD, nIf]);
        % Bootstrap raw correlation between fiber count and delays
        [R_fibNtwD_as(nIf), P_fibNtwD_as(nIf)] = iforest_spearman(fibNas, IEDtwDas);
        % Bootstrap residual correlations between fiber count and delays
        [R_fibNtwD_rED_as(nIf, :), P_fibNtwD_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(fibNas, IEDtwDas, EDas);
    end
    %% Save all results to file
    saveFolder = 'step8_IEDtravelingWave_corrSelect';
    mkdir(saveFolder);
    saveName = ['IEDtw_meth-corrSelect_TD-', TDstr{nTD}, '_polyfit-fibnVSdelay.mat'];
    save(fullfile(saveFolder, saveName), ...
        'R_twDfibN', 'P_twDfibN', 'R_twDfibN_rED', 'P_twDfibN_rED', ...
        'R_fibNtwD_as', 'P_fibNtwD_as', 'R_fibNtwD_rED_as', 'P_fibNtwD_rED_as');
end