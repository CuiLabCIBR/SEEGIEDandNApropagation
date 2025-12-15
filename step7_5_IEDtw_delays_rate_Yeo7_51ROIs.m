% Calculate the relationship between velocity and fiber QA, and between 
% velocity and fiber length.
clc; clear; close all;
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
elecAnatPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDtw_meth-None_TD-6d0.mat'));
%% Main Processing Loop mapping each channel to Yeo 7 network
% Initialize group-level matrices
IEDtwD_Yeo7_51ROIs_as = zeros(51, 51, 47);
IEDtwR_Yeo7_51ROIs_as = zeros(51, 51, 47);
% Process each subject individually
for nS = 1:length(subG)
    subID = subG{nS};
    chanLabel = IEDtw(nS).label;
    %% Extract IED transfer matrices
    IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
    IEDtwR = IEDtw(nS).indirectTransferRateMatrix;
    Mask = IEDtwR;
    Mask = Mask - diag(diag(Mask));
    Mask(isnan(Mask)) = 0;
    Mask(Mask<0.1) = 0;
    IEDtwD(Mask==0) = nan;
    IEDtwR(Mask==0) = nan;
    %% Electrode channel to Brain Region Mapping
    % Load subject-specific electrode anatomy information
    load(fullfile(elecAnatPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']));
    Yeo7Idx = zeros(length(chanLabel), 1);
    % Iterate through all channels
    for nChan = 1:length(chanLabel)
        a = strcmp(chanLabel{nChan}, contAnat.label);
        varList = fieldnames(contAnat);
        if ismember('Yeo7_51_idx', varList)
            if ~isempty(contAnat.Yeo7_51_idx{a==1})
                Yeo7Idx(nChan) = contAnat.Yeo7_51_idx{a==1};
            end
        end
    end
    %% Channel to Region Aggregation
    % Initialize cell arrays for region-level data aggregation
    IEDtwD_Yeo7_51ROIs = cell(51, 51);
    IEDtwR_Yeo7_51ROIs = cell(51, 51);
    % Collect all valid channel pairs between regions
    for nChan1 = 1:size(IEDtwD, 1)
        for nChan2 = 1:size(IEDtwD, 1)
            yeo7Idx1 = Yeo7Idx(nChan1);
            yeo7Idx2 = Yeo7Idx(nChan2);
            % Process only if both channels are mapped to valid regions
            if yeo7Idx1>0 && yeo7Idx2>0
                % Append both connection directions (bidirectional)
                IEDtwD_Yeo7_51ROIs{yeo7Idx1, yeo7Idx2} = [IEDtwD_Yeo7_51ROIs{yeo7Idx1, yeo7Idx2}, IEDtwD(nChan1, nChan2)];
                IEDtwR_Yeo7_51ROIs{yeo7Idx1, yeo7Idx2} = [IEDtwR_Yeo7_51ROIs{yeo7Idx1, yeo7Idx2}, IEDtwR(nChan1, nChan2)];
            end
        end
    end
    %% Subject-Level Region Averages
    % Compute mean values per region pair (ignoring NaNs)
    for n1 = 1:51
        for n2 = 1:51
            IEDtwD_Yeo7_51ROIs_as(n1, n2, nS) = mean(IEDtwD_Yeo7_51ROIs{n1, n2}, 'omitnan');
            IEDtwR_Yeo7_51ROIs_as(n1, n2, nS) = mean(IEDtwR_Yeo7_51ROIs{n1, n2}, 'omitnan');
        end
    end
end
IEDtwD_Yeo7_51ROIs_as_mean = mean(IEDtwD_Yeo7_51ROIs_as, 3, 'omitnan');
IEDtwR_Yeo7_51ROIs_as_mean = mean(IEDtwR_Yeo7_51ROIs_as, 3, 'omitnan');
saveFolder = 'step7_IEDtw_Yeo7';
saveName = 'IEDtw_delays_rate_Yeo7_51ROIs.mat';
save(fullfile(saveFolder, saveName), 'IEDtwD_Yeo7_51ROIs_as_mean', 'IEDtwR_Yeo7_51ROIs_as_mean');