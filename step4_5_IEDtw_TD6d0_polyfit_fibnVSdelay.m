% Analyzes the correlation between fiber count and IED propagation delay 
% while controlling for Euclidean distance effects through polynomial 
% fitting and bootstrap validation.
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
load(fullfile('step2_fiberTrack', 'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDtw_meth-None_TD-6d0.mat'));
for nS = 1:length(subG)
    subID = subG{nS};
    % Load Euclidean distance matrix for current subject
    load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
    %% Extract fiber count data
    fibN = fiberCount{nS};
    if isempty(fibN)
        continue;
    end
    % Extract IED traveling wave data
    IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
    IEDtwR = IEDtw(nS).indirectTransferRateMatrix;
    Mask = IEDtwR;
    Mask = Mask - diag(diag(Mask));
    Mask(isnan(Mask)) = 0;
    Mask(Mask<0.1) = 0;
    Mask(fibN==0) = 0;
    ED = DistanceWorld(Mask>0);
    IEDtwD = IEDtwD(Mask>0);
    fibN = fibN(Mask>0);
    fibN = log(fibN);
    %% Correlation Analysis
    % Raw correlation between fiber count and delay
    [R_fibNtwD(nS), P_fibNtwD(nS)] = corr(fibN, IEDtwD, "type", "Spearman");
    % Remove distance effect from fiber count using polynomial fit
    [polyModel, ~, ~] = f_polyfit(ED, fibN, 6);
    fibN_rED = polyModel.residual;
    % Remove distance effect from delays using polynomial fit
    [polyModel, ~, ~] = f_polyfit(ED, IEDtwD, 6);
    IEDtwD_rED = polyModel.residual;
    % Correlation between residuals (rate vs delay)
    [R_fibNtwD_rED(nS), P_fibNtwD_rED(nS)]  = corr(fibN_rED, IEDtwD_rED, "type", "Spearman");
end
saveFolder = 'step4_IEDtw_polyfit_TD6d0';
mkdir(saveFolder);
saveName = 'IEDtw_meth-none_TD-6d0_polyfit-fibnVSdelay.mat';
save(fullfile(saveFolder, saveName), ...
    'R_fibNtwD', 'P_fibNtwD', 'R_fibNtwD_rED', 'P_fibNtwD_rED');
% %% Aggragate all data 
% EDas = []; IEDtwDas = []; fibNas = [];
% EDas = [EDas; ED]; 
% IEDtwDas = [IEDtwDas; IEDtwD]; 
% fibNas = [fibNas; fibN];
% %% Bootstrap analysis using parallel processing
% parfor nIf = 1:1000
%     disp(nIf);
%     % Bootstrap raw correlation between fiber count and delays
%     [R_fibNtwD_as(nIf), P_fibNtwD_as(nIf)] = iforest_spearman(fibNas, IEDtwDas);
%     % Bootstrap residual correlations between fiber count and delays
%     [R_fibNtwD_rED_as(nIf, :), P_fibNtwD_rED_as(nIf, :)] ...
%         = iforest_spearman_regressConfound(fibNas, IEDtwDas, EDas);
% end
% saveFolder = 'step4_IEDtw_polyfit_TD6d0';
% mkdir(saveFolder);
% saveName = 'IEDtw_meth-none_TD-6d0_polyfit-fibnVSdelay.mat';
% save(fullfile(saveFolder, saveName), ...
%     'EDas', 'IEDtwDas', 'fibNas', ...
%     'R_fibNtwD', 'P_fibNtwD', 'R_fibNtwD_rED', 'P_fibNtwD_rED', ...
%     'R_fibNtwD_as', 'P_fibNtwD_as', 'R_fibNtwD_rED_as', 'P_fibNtwD_rED_as');