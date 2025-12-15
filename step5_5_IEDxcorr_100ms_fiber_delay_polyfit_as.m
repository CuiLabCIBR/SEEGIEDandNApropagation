% Analyzes the relationship between fiber count and IED cross-correlation 
% delays of each subject, while controlling for the effect of 
% inter-channel Euclidean distance using polynomial regression. 
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
load(fullfile('step2_fiberTrack',  'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
% Load IED cross-correlation data
load(fullfile('step5_IEDxcorr_polyfit', 'IED6d0_xcorr100ms_Rmax_delay_allSubject.mat'), ...
    'Delay_allSubject', 'xcR_allSubject');
%% Process each subject
% Initialize aggregated data storage
IEDxcDas = []; EDas = []; fibNas = [];
for nS = 1:length(subG)
        subID = subG{nS};
        % Load inter-channel Euclidean distance matrix
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']), 'DistanceWorld');
        % Extract subject-specific data
        IEDxcD = Delay_allSubject{nS};
        IEDxcR = xcR_allSubject{nS};
        fibN = fiberCount{nS};
        if isempty(fibN); continue; end
        % Extract valid data
        Mask = IEDxcR;
        Mask = Mask - diag(diag(Mask));
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(IEDxcR) | isinf(IEDxcR) | IEDxcR <= zThd) = 0;
        Mask(isnan(IEDxcD) | isinf(IEDxcD)) = 0;
        Mask(fibN==0) = 0;
        ED = DistanceWorld(Mask>0);
        IEDxcD = IEDxcD(Mask>0);
        fibN = log(fibN(Mask>0));
       % Aggregate data across subjects
        EDas = [EDas; ED];
        IEDxcDas = [IEDxcDas; IEDxcD];
        fibNas = [fibNas; fibN];
end
 % Compute correlation between fiber count and IED cross-correlation delays
[R_fibNxcD_as, P_fibNxcD_as] = corr(fibNas, IEDxcDas, "type", "Spearman");
% Compute distance-residualized correlations
% Remove distance effect from delays using polynomial fit
[polyModel, ~, ~] = f_polyfit(EDas, IEDxcDas, 6);
IEDxcDres = polyModel.residual+polyModel.p(end);
% Remove distance effect from fiber counts using polynomial fit
[polyModel, ~, ~] = f_polyfit(EDas, fibNas, 6);
fibNres = polyModel.residual+polyModel.p(end);
% Compute correlation between residuals
[R_fibNxcD_rED_as, P_fibNxcD_rED_as] = corr(fibNres, IEDxcDres, "type", "Spearman");
% Save results
saveFolder = 'step5_IEDxcorr_polyfit';
saveName = 'IED6d0_XC100ms_polyfit_delay_fibN_as';
save(fullfile(saveFolder, saveName), ...
    'fibNas', 'IEDxcDas',  'IEDxcDres', 'fibNres', ...
    'R_fibNxcD_as', 'P_fibNxcD_as', ...
    'R_fibNxcD_rED_as', 'P_fibNxcD_rED_as');