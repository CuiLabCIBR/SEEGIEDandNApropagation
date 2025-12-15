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
%% Load data for each subject and each null model iteration
% This section loads IED traveling wave null model data and Euclidean distance
% matrices for all subjects, organizing them for subsequent analysis
for nS = 1:length(subG)
    subID = subG{nS};
    % Load IED traveling wave null model data for current subject
    load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDtwNull', [subID, '_IEDtwNull_meth-None_TD-6d0.mat']));
    % Load Euclidean distance matrix for current subject
    load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
    distanceAll{nS} = DistanceWorld;
    for nNull = 1:1000
        delayAll{nNull}{nS} = IEDtwNull.indirectTransferDelayMatrix(:, :, nNull);
        rateAll{nNull}{nS} = IEDtwNull.indirectTransferRateMatrix(:, :, nNull);
    end
end
%% Polynomial fitting analysis for null models
% Initialize arrays to store correlation results and fit quality metrics
R_twDR = zeros(length(subG), 1000); 
P_twDR = zeros(length(subG), 1000);
R_twDR_rED = zeros(length(subG), 6, 1000); 
P_twDR_rED = zeros(length(subG), 6, 1000);
adjR2_twDED = zeros(length(subG), 6, 1000); 
parfor nNull = 1:1000
    disp(nNull);
    delayNull = delayAll{nNull};
    rateNull = rateAll{nNull};
    op  = subparfor(delayNull, rateNull, distanceAll);
    % subject level polyfit
    R_twDR(:, nNull) = op.R_twDR; 
    P_twDR(:, nNull) = op.P_twDR;
    R_twDR_rED(:, :, nNull) = op.R_twDR_rED; 
    P_twDR_rED(:, :, nNull) = op.P_twDR_rED;
    adjR2_twDED(:, :, nNull) = op.adjR2_twDED;
end
saveFolder = 'step4_IEDtw_polyfit_TD6d0';
mkdir(saveFolder);
saveName = 'IEDtwNull_meth-none_TD-6d0_polyfit_distanceVSdelay_rateVSdelay.mat';
save(fullfile(saveFolder, saveName),  'adjR2_twDED', ...
    'R_twDR', 'P_twDR', 'R_twDR_rED', 'P_twDR_rED');
%% subfunction for parfor
function op = subparfor(delayNull, rateNull, distanceAll)
    for nS = 1:47
        IEDtwD = delayNull{nS};
        IEDtwR = rateNull{nS};
        DistanceWorld = distanceAll{nS};
        Mask = IEDtwR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask<0.1) = 0;
        ED = DistanceWorld(Mask>0);
        IEDtwD = IEDtwD(Mask>0);
        IEDtwR = IEDtwR(Mask>0);
        IEDtwR = log(IEDtwR);
        % Correlation Analysis
        [R_twDR(nS), P_twDR(nS)] = corr(IEDtwR, IEDtwD, "type", "Spearman");
        for order = 1:6
            [polyModel, ~, ~] = f_polyfit(ED, IEDtwR, order);
            IEDtwR_rED = polyModel.residual;
            [polyModel, ~, ~] = f_polyfit(ED, IEDtwD, order);
            IEDtwD_rED = polyModel.residual;
            adjR2_twDED(nS, order) = polyModel.adjR2;
            [R_twDR_rED(nS, order), P_twDR_rED(nS, order)] = corr(IEDtwR_rED, IEDtwD_rED, "type", "Spearman");
        end
    end
    op.R_twDR = R_twDR; 
    op.P_twDR = P_twDR;
    op.R_twDR_rED = R_twDR_rED; 
    op.P_twDR_rED = P_twDR_rED;
    op.adjR2_twDED = adjR2_twDED;
end
%% aggregate analysis
% IfN = 25;
% R_twDR_as = zeros(IfN, 1000); 
% P_twDR_as = zeros(IfN, 1000);
% R_twDR_rED_as = zeros(IfN, 6, 1000); 
% P_twDR_rED_as = zeros(IfN, 6, 1000);
% adjR2_twDED_as = zeros(IfN, 6, 1000);
%         EDas = []; 
%     IEDtwDas = []; 
%     IEDtwRas = [];
%      % Aggragate all data 
%         EDas = [EDas; ED]; IEDtwDas = [IEDtwDas; IEDtwD]; IEDtwRas = [IEDtwRas; IEDtwR];
%     for nIf = 1:25
%         % Bootstrap raw correlation between delays and rate
%         [R_twDR_as(nIf), P_twDR_as(nIf)] = iforest_spearman(IEDtwDas, IEDtwRas);
%         % Bootstrap residual correlations between delays and rate
%         [R_twDR_rED_as(nIf, :), P_twDR_rED_as(nIf, :)] ...
%             = iforest_spearman_regressConfound(IEDtwDas, IEDtwRas, EDas);
%         % Bootstrap polynomial fit quality between Euclidean distance and
%         % delays
%         adjR2_twDED_as(nIf, :) = iforest_polyfit(EDas, IEDtwDas);
%     end
%     op2.R_twDR_as = R_twDR_as; 
%     op2.P_twDR_as = P_twDR_as;
%     op2.R_twDR_rED_as = R_twDR_rED_as; 
%     op2.P_twDR_rED_as = P_twDR_rED_as;
%     op2.adjR2_twDED_as = adjR2_twDED_as;