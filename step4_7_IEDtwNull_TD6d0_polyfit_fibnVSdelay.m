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
load(fullfile('step2_fiberTrack', 'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
% Initialize cell arrays to store delay, rate, and distance data
for nS = 1:length(subG)
    subID = subG{nS};
    % Load IEDtwNull data for each subject
    load(fullfile('step4_IEDtw_polyfit_TD6d0', ...
        'IEDtwNull', [subID, '_IEDtwNull_meth-None_TD-6d0.mat']));
    for nNull = 1:1000
        delayAll{nNull}{nS} = IEDtwNull.indirectTransferDelayMatrix(:, :, nNull);
        rateAll{nNull}{nS} = IEDtwNull.indirectTransferRateMatrix(:, :, nNull);
    end
    % Load Euclidean distance matrix for each subject
    load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
    fibN{nS} = fiberCount{nS};
    EDall{nS} = DistanceWorld;
end
% Initialize matrices to store correlation results
R_fibNtwD = zeros(length(subG), 1000); 
P_fibNtwD = zeros(length(subG), 1000);
R_fibNtwD_rED = zeros(length(subG), 6, 1000); 
P_fibNtwD_rED = zeros(length(subG), 6, 1000);
IfN = 25;
R_fibNtwD_as = zeros(IfN, 1000); 
P_fibNtwD_as = zeros(IfN, 1000);
R_fibNtwD_rED_as = zeros(IfN, 6, 1000); 
P_fibNtwD_rED_as = zeros(IfN, 6, 1000);
% Perform parallel computation for 1000 null iterations
parfor nNull = 1:1000
    disp(nNull);
    delayNull = delayAll{nNull};
    rateNull = rateAll{nNull};
    [R_fibNtwD(:, nNull), P_fibNtwD(:, nNull), ...
        R_fibNtwD_rED(:, :, nNull), P_fibNtwD_rED(:, :, nNull), ...
        R_fibNtwD_as(:, nNull), P_fibNtwD_as(:, nNull), ...
        R_fibNtwD_rED_as(:, :, nNull), P_fibNtwD_rED_as(:, :, nNull)] ...
    = subparfor(delayNull, rateNull, EDall, fibN);
end
% Save results
saveFolder = 'step4_IEDtw_polyfit_TD6d0';
mkdir(saveFolder);
saveName = 'IEDtwNull_meth-none_TD-6d0_polyfit-fibnVSdelay.mat';
save(fullfile(saveFolder, saveName), ...
    'R_fibNtwD', 'P_fibNtwD', 'R_fibNtwD_rED', 'P_fibNtwD_rED', ...
    'R_fibNtwD_as', 'P_fibNtwD_as', 'R_fibNtwD_rED_as', 'P_fibNtwD_rED_as');
%% Subfunction for parallel processing
function [R_fibNtwD, P_fibNtwD, R_fibNtwD_rED, P_fibNtwD_rED, ...
    R_fibNtwD_as, P_fibNtwD_as, R_fibNtwD_rED_as, P_fibNtwD_rED_as] ...
    = subparfor(IEDtwD_null, IEDtwR_null, EDall, fibN)
    EDas = []; IEDtwDas = []; fibNas = [];
    R_fibNtwD = zeros(length(IEDtwD_null), 1); 
    P_fibNtwD = zeros(length(IEDtwD_null), 1);
    R_fibNtwD_rED = zeros(length(IEDtwD_null), 6); 
    P_fibNtwD_rED = zeros(length(IEDtwD_null), 6);
    % Process each subject's data
    for nS = 1:length(IEDtwD_null)
        fibN1 = fibN{nS};
        if ~isempty(fibN1)
            IEDtwD = IEDtwD_null{nS};
            IEDtwR = IEDtwR_null{nS};
            ED = EDall{nS};
            % Apply mask to extract relevant data points
            Mask = IEDtwR;
            Mask = Mask - diag(diag(Mask));
            Mask(isnan(Mask)) = 0;
            Mask(Mask<0.1) = 0;
            Mask(fibN1==0) = 0;
            ED = ED(Mask>0);
            IEDtwD = IEDtwD(Mask>0);
            fibN1 = fibN1(Mask>0);
            fibN1 = log(fibN1);
            % Correlation Analysis: fiber count vs. delay
            [R_fibNtwD(nS), P_fibNtwD(nS)] = corr(fibN1, IEDtwD, "type", "Spearman");
            for order = 1:6
                % Residual correlation after removing Euclidean distance effect
                [polyModel, ~, ~] = f_polyfit(ED, fibN1, order);
                fibN_rED = polyModel.residual;
                [polyModel, ~, ~] = f_polyfit(ED, IEDtwD, order);
                IEDtwD_rED = polyModel.residual;
                [R_fibNtwD_rED(nS, order), P_fibNtwD_rED(nS, order)] = corr(fibN_rED, IEDtwD_rED, "type", "Spearman");
            end
            % Aggragate all data 
            EDas = [EDas; ED]; 
            IEDtwDas = [IEDtwDas; IEDtwD]; 
            fibNas = [fibNas; fibN1];
        end
    end
    % Initialize arrays for bootstrap analysis
    R_fibNtwD_as = zeros(25, 1); 
    P_fibNtwD_as = zeros(25, 1);
    R_fibNtwD_rED_as = zeros(25, 6); 
    P_fibNtwD_rED_as = zeros(25, 6);
    for nIf = 1:25
        % Bootstrap raw correlation between fiber count and delays
        [R_fibNtwD_as(nIf), P_fibNtwD_as(nIf)] = iforest_spearman(fibNas, IEDtwDas);
        % Bootstrap residual correlations between fiber count and delays
        [R_fibNtwD_rED_as(nIf, :), P_fibNtwD_rED_as(nIf, :)] ...
            = iforest_spearman_regressConfound(fibNas, IEDtwDas, EDas);
    end
end