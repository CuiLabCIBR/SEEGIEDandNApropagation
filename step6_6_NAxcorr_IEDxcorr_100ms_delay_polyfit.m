% Analyzes the relationship between IED and non-IED cross-correlation delays
% across multiple subjects using Spearman correlation while controlling 
% for Euclidean distance effects through polynomial regression (orders 1-6). 
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
% Load IED and NA cross-correlation data
load(fullfile('step6_NAxcorr_polyfit',  'NAnon6d0IED_xcorr100ms_Rmax_delay_allSubject.mat'));
NAdelay = Delay_allSubject;
NAr = xcR_allSubject;
load(fullfile('step5_IEDxcorr_polyfit', 'IED6d0_xcorr100ms_Rmax_delay_allSubject.mat'));
IEDdelay = Delay_allSubject;
IEDr = xcR_allSubject;
EDas = []; IEDxcDas = []; NAxcDas = [];
% Process each subject individually
for nS = 1:length(subG)
    subID = subG{nS};
    % Load Euclidean distance matrix
    load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']));
    %% Extract transfer matrix
    % Extract IED cross-correlation delay and correlation coefficient matrices
    IEDxcD = IEDdelay{nS};
    IEDxcR = IEDr{nS};
    % Extract NA cross-correlation delay and correlation coefficient matrices
    NAxcD = NAdelay{nS};
    NAxcR = NAr{nS};
    % Extract valid data points using the mask
    Mask = IEDxcR;
    Mask = Mask - diag(diag(Mask));
    zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
    Mask(isnan(NAxcR) | isinf(NAxcR) | NAxcR <= zThd) = 0;
    Mask(isnan(NAxcD) | isinf(NAxcD)) = 0;
    Mask(isnan(IEDxcR) | isinf(IEDxcR) | IEDxcR <= zThd) = 0;
    Mask(isnan(IEDxcD) | isinf(IEDxcD)) = 0;
    ED = DistanceWorld(Mask>0);
    IEDxcD = IEDxcD(Mask>0);
    NAxcD = NAxcD(Mask>0);
    EDas = [EDas; ED];
    IEDxcDas = [IEDxcDas; IEDxcD];
    NAxcDas = [NAxcDas; NAxcD];
    %% Correlation analysis
    [RxcD(nS), PxcD(nS)] = corr(IEDxcD, NAxcD, 'type', 'Spearman');
    for order = 1:6
        [polyModel, ~, ~] = f_polyfit(ED, NAxcD, order);
        NAxcDres = polyModel.residual;
        [polyModel, ~, ~] = f_polyfit(ED, IEDxcD, order);
        IEDxcDres = polyModel.residual;
        [RxcD_rED(nS, order), PxcD_rED(nS, order)] ...
            = corr(NAxcDres, IEDxcDres, 'type', 'Spearman');
    end
end
 %% Parallel processing for bootstrap analysis (1000 iterations)
parfor nIf = 1:1000
    disp(nIf);
    [RxcD_as(nIf), PxcD_as(nIf)] = iforest_spearman(IEDxcDas, NAxcDas);
    [RxcD_rED_as(nIf, :), PxcD_rED_as(nIf, :)] ...
        = iforest_spearman_regressConfound(IEDxcDas, NAxcDas, EDas);
end
%% Save all analysis results
saveFolder = 'step6_NAxcorr_polyfit';
saveName = 'NAnon6d0IED_XC100ms_IEDxcD_polyfit_delays.mat';
save(fullfile(saveFolder, saveName), ...
    "RxcD", "PxcD",  "RxcD_rED", "PxcD_rED", ...
    "EDas", "IEDxcDas", "NAxcDas", ...
     "RxcD_as", "PxcD_as", "RxcD_rED_as", "PxcD_rED_as");